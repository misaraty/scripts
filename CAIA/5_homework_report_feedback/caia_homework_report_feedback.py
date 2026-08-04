import os
import sys

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import json
from pathlib import Path

import matplotlib.pyplot as plt

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

HOMEWORK_TEXT = """
乙酸乙酯皂化是二级反应。由速率方程积分后可以用电导率随时间的变化求速率常数。
我认为拟合曲线的R2很高，所以模型一定正确。活化能用两个温度下的速率常数代入阿伦尼乌斯公式计算。
""".strip()

REPORT_TEXT = """
实验目的：测定乙酸乙酯皂化反应的速率常数和活化能。
实验原理：反应过程中OH-被CH3COO-取代，电导率逐渐降低。
数据处理：使用非线性拟合得到28℃下k=8.58，38℃下k=18.37，并计算活化能。
结果分析：两条曲线拟合较好。误差可能来自混合不充分、恒温控制和电导率读数。
结论：获得了两个温度下的速率常数。
""".strip()

HOMEWORK_RULES = [
    ("反应级数", ["二级反应"], 20),
    ("积分速率方程", ["积分", "速率方程"], 20),
    ("实验观测量", ["电导率"], 15),
    ("模型评价", ["残差", "RMSE", "MAE", "物理意义"], 25),
    ("活化能", ["阿伦尼乌斯", "两个温度"], 20),
]

REPORT_RULES = [
    ("实验目的", ["实验目的", "速率常数", "活化能"], 12.5),
    ("实验原理", ["实验原理", "OH", "电导率"], 12.5),
    ("原始数据与单位", ["min", "μS", "mol", "单位"], 12.5),
    ("数据处理", ["拟合", "k="], 12.5),
    ("参数解释", ["参数", "速率常数", "物理意义"], 12.5),
    ("残差分析", ["残差", "RMSE", "最大绝对残差"], 12.5),
    ("误差分析", ["误差", "恒温", "混合", "读数"], 12.5),
    ("结论表达", ["结论", "两个温度", "活化能"], 12.5),
]


def score_text(text: str, rules: list[tuple[str, list[str], float]]) -> list[dict]:
    lower = text.lower()
    results = []
    for dimension, keywords, full_score in rules:
        hits = [keyword for keyword in keywords if keyword.lower() in lower]
        ratio = len(hits) / len(keywords)
        score = full_score * min(1.0, ratio * 1.25)
        results.append(
            {
                "dimension": dimension,
                "score": round(score, 2),
                "full_score": full_score,
                "hit_keywords": hits,
                "missing_keywords": [
                    keyword for keyword in keywords if keyword not in hits
                ],
            }
        )
    return results


def homework_feedback(scores: list[dict]) -> str:
    total = sum(item["score"] for item in scores)
    return (
        f"理论作业演示评分为{total:.1f}/100。该分数只反映预设关键词的覆盖情况，不代表程序已经判断了全部内容是否正确。\n\n"
        "作答说明了乙酸乙酯皂化的反应级数、电导率观测量以及阿伦尼乌斯公式的用途。关于模型评价，R²较高只能说明拟合曲线对当前数据的解释程度较高，不能据此认定模型一定正确。还需要结合残差、原始尺度上的RMSE或MAE、端点参数及其物理意义进行判断。\n\n"
        "可继续比较28.0 ℃下三种处理方法得到的速率常数，并说明两温度法计算活化能时为何会受到速率常数比值的影响。"
    )


def report_feedback(scores: list[dict]) -> str:
    total = sum(item["score"] for item in scores)
    weak = sorted(scores, key=lambda item: item["score"] / item["full_score"])[:3]
    weak_names = "、".join(item["dimension"] for item in weak)
    missing_terms = []
    for item in weak:
        missing_terms.extend(item["missing_keywords"])
    missing_text = "、".join(dict.fromkeys(missing_terms))

    return (
        f"实验报告演示评分为{total:.1f}/100。该分数来自关键词覆盖规则，仅用于展示形成性反馈流程。\n\n"
        "报告已经写明实验目的、反应原理、两个温度下的拟合结果以及主要误差来源。"
        f"目前较薄弱的部分是{weak_names}。可补充原始数据及单位、参数的物理意义和残差评价，并在结论中同时给出两个速率常数、活化能和主要误差来源。"
        + (f"当前文本中未检出这些示例关键词：{missing_text}。" if missing_text else "")
    )


def optional_deepseek_polish(base_feedback: str, original_text: str) -> str | None:
    enabled = os.getenv("CAIA_USE_DEEPSEEK", "0").strip() == "1"
    api_key = os.getenv("DEEPSEEK_API_KEY", "XXXXXX").strip()
    model = os.getenv("DEEPSEEK_MODEL", "deepseek-v4-pro").strip()
    base_url = os.getenv("DEEPSEEK_BASE_URL", "https://api.deepseek.com").rstrip("/")
    if not enabled or not api_key or not model:
        return None

    first_paragraph, separator, remainder = base_feedback.partition("\n\n")
    body = remainder if separator else base_feedback

    try:
        import requests

        prompt = (
            "请把下面的反馈正文改写得清楚、自然、简洁。不得修改评价结论，不得新增事实或数值。"
            "使用连续自然段，不使用模板化小标题、项目符号、引号或破折号。不要重复评分。\n\n"
            f"学生文本：\n{original_text}\n\n反馈正文：\n{body}"
        )
        response = requests.post(
            f"{base_url}/chat/completions",
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json",
            },
            json={
                "model": model,
                "messages": [{"role": "user", "content": prompt}],
                "temperature": 0.2,
                "stream": False,
            },
            timeout=120,
        )
        response.raise_for_status()
        polished_body = response.json()["choices"][0]["message"]["content"].strip()
        return f"{first_paragraph}\n\n{polished_body}" if polished_body else None
    except Exception as exc:
        print(f"DeepSeek调用失败：{type(exc).__name__}", file=sys.stderr)
        return None


def draw_dimension_chart(report_scores: list[dict]) -> None:
    labels = [item["dimension"] for item in report_scores]
    percentages = [100 * item["score"] / item["full_score"] for item in report_scores]
    plt.figure(figsize=(9.2, 5.4))
    bars = plt.bar(labels, percentages)
    for bar, value in zip(bars, percentages):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value + 1.5,
            f"{value:.0f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    plt.ylabel("维度达成率 / %", fontsize=11)
    plt.title("实验报告关键词覆盖情况", fontsize=14)
    plt.ylim(0, 110)
    plt.xticks(rotation=25, ha="right")
    plt.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "实验报告评价维度.jpg", dpi=600, format="jpg", bbox_inches="tight"
    )
    plt.close()


def main() -> None:
    homework_scores = score_text(HOMEWORK_TEXT, HOMEWORK_RULES)
    report_scores = score_text(REPORT_TEXT, REPORT_RULES)

    base_homework = homework_feedback(homework_scores)
    base_report = report_feedback(report_scores)
    homework_text = (
        optional_deepseek_polish(base_homework, HOMEWORK_TEXT) or base_homework
    )
    report_text = optional_deepseek_polish(base_report, REPORT_TEXT) or base_report
    feedback_mode = (
        "deepseek_polish"
        if homework_text != base_homework or report_text != base_report
        else "rule_based"
    )

    payload = {
        "homework": {
            "original": HOMEWORK_TEXT,
            "scores": homework_scores,
            "base_feedback": base_homework,
            "feedback": homework_text,
        },
        "lab_report": {
            "original": REPORT_TEXT,
            "scores": report_scores,
            "base_feedback": base_report,
            "feedback": report_text,
        },
        "feedback_mode": feedback_mode,
        "notice": "评分依据是预设关键词的覆盖比例。",
    }
    (OUTPUT_DIR / "作业与实验报告分析.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (OUTPUT_DIR / "理论作业反馈.txt").write_text(homework_text, encoding="utf-8")
    (OUTPUT_DIR / "实验报告反馈.txt").write_text(report_text, encoding="utf-8")
    draw_dimension_chart(report_scores)
    print(homework_text)
    print("\n" + report_text)
    print("\n输出目录：", OUTPUT_DIR.resolve())


if __name__ == "__main__":
    main()
