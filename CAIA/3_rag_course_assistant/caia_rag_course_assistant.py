import os
import sys

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

try:
    from dotenv import load_dotenv

    load_dotenv(Path(__file__).resolve().parents[1] / ".env")
    load_dotenv()
except ImportError:
    pass

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Microsoft YaHei",
    "Arial",
    "Noto Sans CJK JP",
    "Noto Sans CJK SC",
    "SimHei",
    "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False

OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

# 公开仓库只放原创演示片段。正式使用时应替换为获得授权的教材、讲义、习题和实验规范。
DOCUMENTS = [
    {
        "source": "动力学理论讲义",
        "section": "一级反应",
        "text": "一级反应满足-dc/dt=kc，积分后得到ln(c/c0)=-kt。一级反应半衰期t1/2=ln2/k，与初始浓度无关。判断反应级数时必须同时考察积分式、速率常数单位和实验条件。",
    },
    {
        "source": "动力学理论讲义",
        "section": "二级反应",
        "text": "等初始浓度二级反应满足dx/dt=k(a-x)^2，积分后得到x/[a(a-x)]=kt。二级反应速率常数单位通常为L·mol^-1·time^-1，半衰期与初始浓度有关。",
    },
    {
        "source": "物理化学实验讲义",
        "section": "蔗糖水解",
        "text": "蔗糖水解在水和氢离子浓度近似恒定时可视为准一级反应。旋光度与体系组成相关，可写成alpha_t=alpha_inf+(alpha_0-alpha_inf)exp(-kt)，由此求得速率常数和半衰期。",
    },
    {
        "source": "物理化学实验讲义",
        "section": "乙酸乙酯皂化",
        "text": "乙酸乙酯与氢氧化钠等浓度反应可按二级反应处理。反应中高电导率的OH-逐渐被低电导率的CH3COO-取代，电导率随时间降低，模型可写为kappa_t=B/(1+At)+C，其中A=c0k。",
    },
    {
        "source": "数据处理规范",
        "section": "模型评价",
        "text": "模型评价不能只看R2。应在原始实验信号尺度上综合比较RMSE、MAE、最大绝对残差、参数稳定性和物理意义。曲线在图像上接近，并不表示参数完全一致。",
    },
    {
        "source": "实验报告评价量规",
        "section": "反馈原则",
        "text": "实验报告反馈应检查实验原理、原始数据和单位、数据处理、参数物理意义、残差与误差分析以及结论表达。人工智能反馈用于形成性评价，最终判断由教师审核。",
    },
    {
        "source": "课程人工智能使用规范",
        "section": "可信使用",
        "text": "确定性程序负责数值计算，RAG负责检索课程证据，大语言模型负责语言解释和练习生成，教师负责审核。证据不足时应明确说明，不得编造教材页码、数据或引用。",
    },
]

DEFAULT_QUESTION = "一级反应半衰期与初始浓度有什么关系？"
REACTION_ORDERS = ("零级反应", "一级反应", "二级反应")


def tokenize(text: str) -> list[str]:
    text = text.lower()
    latin = re.findall(r"[a-z]+\d*|\d+(?:\.\d+)?", text)
    chinese_tokens: list[str] = []
    for segment in re.findall(r"[\u4e00-\u9fff]+", text):
        chinese_tokens.extend(segment)
        chinese_tokens.extend(segment[i : i + 2] for i in range(len(segment) - 1))
    return latin + chinese_tokens


class BM25:
    def __init__(self, documents: list[str], k1: float = 1.5, b: float = 0.75):
        self.tokens = [tokenize(doc) for doc in documents]
        self.k1 = k1
        self.b = b
        self.lengths = np.array([len(x) for x in self.tokens], dtype=float)
        self.avgdl = float(np.mean(self.lengths)) if len(self.lengths) else 1.0
        self.term_freq = [Counter(x) for x in self.tokens]
        document_frequency = Counter()
        for tokens in self.tokens:
            document_frequency.update(set(tokens))
        n = len(self.tokens)
        self.idf = {
            term: math.log(1.0 + (n - frequency + 0.5) / (frequency + 0.5))
            for term, frequency in document_frequency.items()
        }

    def score(self, query: str) -> np.ndarray:
        query_tokens = tokenize(query)
        scores = np.zeros(len(self.tokens), dtype=float)
        for i, tf in enumerate(self.term_freq):
            dl = self.lengths[i]
            for term in query_tokens:
                if term not in tf:
                    continue
                frequency = tf[term]
                denominator = frequency + self.k1 * (
                    1.0 - self.b + self.b * dl / max(self.avgdl, 1e-9)
                )
                scores[i] += (
                    self.idf.get(term, 0.0) * frequency * (self.k1 + 1.0) / denominator
                )
        return scores


def minmax(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0 or np.allclose(values.max(), values.min()):
        return np.zeros_like(values)
    return (values - values.min()) / (values.max() - values.min())


def detect_reaction_order(text: str) -> str | None:
    for order in REACTION_ORDERS:
        if order in text:
            return order
    return None


def topic_adjustment(question: str, section: str, text: str) -> float:
    target_order = detect_reaction_order(question)
    document_text = f"{section} {text}"
    adjustment = 0.0

    if section and section in question:
        adjustment += 0.20

    if target_order:
        if target_order in document_text:
            adjustment += 0.35
        elif any(
            order in document_text for order in REACTION_ORDERS if order != target_order
        ):
            adjustment -= 0.35

    return adjustment


def retrieve(question: str, top_k: int = 4, vector_weight: float = 0.55) -> list[dict]:
    corpus = [f"{item['section']} {item['text']}" for item in DOCUMENTS]
    vectorizer = TfidfVectorizer(analyzer="char", ngram_range=(2, 4), min_df=1)
    matrix = vectorizer.fit_transform(corpus + [question])
    vector_scores = cosine_similarity(matrix[-1], matrix[:-1]).ravel()
    bm25_scores = BM25(corpus).score(question)
    base_scores = vector_weight * minmax(vector_scores) + (
        1.0 - vector_weight
    ) * minmax(bm25_scores)
    adjustments = np.array(
        [
            topic_adjustment(question, item["section"], item["text"])
            for item in DOCUMENTS
        ],
        dtype=float,
    )
    combined = np.clip(base_scores + adjustments, 0.0, 1.0)
    order = np.argsort(combined)[::-1][: min(max(1, top_k), len(DOCUMENTS))]

    results = []
    for rank, index in enumerate(order, start=1):
        item = DOCUMENTS[int(index)]
        results.append(
            {
                "rank": rank,
                "source": item["source"],
                "section": item["section"],
                "text": item["text"],
                "vector_score": float(vector_scores[index]),
                "bm25_score": float(bm25_scores[index]),
                "base_score": float(base_scores[index]),
                "topic_adjustment": float(adjustments[index]),
                "combined_score": float(combined[index]),
            }
        )
    return results


def select_answer_evidence(
    question: str, evidence: list[dict], max_items: int = 2
) -> list[dict]:
    target_order = detect_reaction_order(question)
    selected: list[dict] = []
    item_limit = 1 if target_order else max_items

    for item in evidence:
        if item["combined_score"] < 0.10:
            continue
        if target_order:
            document_text = f"{item['section']} {item['text']}"
            contains_other_order = any(
                order in document_text
                for order in REACTION_ORDERS
                if order != target_order
            )
            if contains_other_order and target_order not in document_text:
                continue
        selected.append(item)
        if len(selected) >= item_limit:
            break

    return selected


def offline_answer(question: str, evidence: list[dict]) -> tuple[str, list[dict]]:
    selected = select_answer_evidence(question, evidence)
    if not selected:
        return "当前演示知识库证据不足，建议补充课程资料或由教师进一步核查。", []

    paragraphs = [f"问题：{question}", ""]
    for item in selected:
        paragraphs.append(
            f"{item['text']}（来源{item['rank']}：{item['source']}，{item['section']}）"
        )
    paragraphs.extend(
        [
            "",
            "以上内容由课程知识库检索结果整理。涉及具体实验条件时，还需结合课程讲义和实验规范。",
        ]
    )
    return "\n".join(paragraphs), selected


def deepseek_answer(question: str, evidence: list[dict]) -> str | None:
    enabled = os.getenv("CAIA_USE_DEEPSEEK", "0").strip() == "1"
    api_key = os.getenv("DEEPSEEK_API_KEY", "XXXXXX").strip()
    model = os.getenv("DEEPSEEK_MODEL", "deepseek-v4-pro").strip()
    base_url = os.getenv("DEEPSEEK_BASE_URL", "https://api.deepseek.com").rstrip("/")
    if not enabled or not api_key or not model or not evidence:
        return None

    try:
        import requests

        evidence_text = "\n\n".join(
            f"[来源{item['rank']}] {item['source']}，{item['section']}\n{item['text']}"
            for item in evidence
        )
        prompt = (
            "仅依据下列课程证据回答问题。不得编造来源、页码或数据；证据不足时明确说明。"
            "使用自然段表达，不使用模板化小标题，不使用项目符号。结尾列出引用来源编号。\n\n"
            f"问题：{question}\n\n证据：\n{evidence_text}"
        )
        response = requests.post(
            f"{base_url}/chat/completions",
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json",
            },
            json={
                "model": model,
                "messages": [
                    {
                        "role": "system",
                        "content": "你是物理化学课程助教。回答必须以给定证据为依据，并保持简洁。",
                    },
                    {"role": "user", "content": prompt},
                ],
                "temperature": 0.2,
                "stream": False,
            },
            timeout=120,
        )
        response.raise_for_status()
        return response.json()["choices"][0]["message"]["content"].strip()
    except Exception as exc:
        print(f"DeepSeek调用失败：{type(exc).__name__}", file=sys.stderr)
        return None


def save_outputs(
    question: str,
    evidence: list[dict],
    used_evidence: list[dict],
    answer: str,
    answer_mode: str,
) -> None:
    payload = {
        "question": question,
        "retrieved_evidence": evidence,
        "answer_evidence": used_evidence,
        "answer_mode": answer_mode,
        "answer": answer,
    }
    (OUTPUT_DIR / "RAG检索与回答.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (OUTPUT_DIR / "RAG回答.txt").write_text(answer, encoding="utf-8")

    labels = [f"来源{x['rank']}\n{x['section']}" for x in evidence]
    values = [x["combined_score"] for x in evidence]
    plt.figure(figsize=(8.0, 5.0))
    bars = plt.bar(labels, values)
    for bar, value in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.015,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    plt.ylabel("混合检索得分", fontsize=11)
    plt.title("RAG混合检索结果", fontsize=14)
    plt.ylim(0, max(values + [0.1]) * 1.20)
    plt.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "RAG检索得分.jpg", dpi=600, format="jpg", bbox_inches="tight"
    )
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="CAIA物理化学课程RAG知识库")
    parser.add_argument(
        "--question", default=DEFAULT_QUESTION, help="需要检索和回答的问题"
    )
    parser.add_argument("--top-k", type=int, default=4, help="返回证据数量")
    args = parser.parse_args()

    evidence = retrieve(args.question, top_k=max(1, args.top_k))
    answer, used_evidence = offline_answer(args.question, evidence)
    answer_mode = "offline"

    llm_answer = deepseek_answer(args.question, used_evidence)
    if llm_answer:
        answer = llm_answer
        answer_mode = "deepseek"

    save_outputs(args.question, evidence, used_evidence, answer, answer_mode)
    print(answer)
    print("\n输出目录：", OUTPUT_DIR.resolve())


if __name__ == "__main__":
    main()
