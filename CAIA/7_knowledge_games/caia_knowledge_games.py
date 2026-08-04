import os

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

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

QUESTIONS = [
    {
        "category": "反应级数判断",
        "question": "一级反应半衰期与初始浓度的关系是？",
        "options": ["成正比", "成反比", "无关", "无法判断"],
        "answer": 2,
        "explanation": "一级反应t1/2=ln2/k，与初始浓度无关。",
    },
    {
        "category": "公式配对",
        "question": "等初始浓度二级反应的积分关系是？",
        "options": ["c=c0-kt", "ln(c/c0)=-kt", "x/[a(a-x)]=kt", "c=c0exp(kt)"],
        "answer": 2,
        "explanation": "等初始浓度二级反应积分后得到x/[a(a-x)]=kt。",
    },
    {
        "category": "实验观测",
        "question": "蔗糖水解实验中用于跟踪反应进程的主要观测量是？",
        "options": ["旋光度", "电导率", "吸光度", "压力"],
        "answer": 0,
        "explanation": "蔗糖及其水解产物旋光性不同，因此可用旋光度跟踪反应。",
    },
    {
        "category": "实验观测",
        "question": "乙酸乙酯皂化过程中电导率降低的主要原因是？",
        "options": ["温度降低", "OH-被CH3COO-取代", "乙醇导电", "水分蒸发"],
        "answer": 1,
        "explanation": "高电导率OH-逐渐被较低电导率CH3COO-取代。",
    },
    {
        "category": "模型评价",
        "question": "拟合曲线看起来重合，能否说明模型参数完全相同？",
        "options": ["一定能", "不能，还需比较参数和残差", "只看R2即可", "只看斜率即可"],
        "answer": 1,
        "explanation": "图像接近不等于参数等价，应比较残差、参数和物理意义。",
    },
    {
        "category": "活化能",
        "question": "两温度法计算活化能时，结果对什么量较敏感？",
        "options": ["两个速率常数之和", "两个速率常数之比", "初始体积", "图像颜色"],
        "answer": 1,
        "explanation": "阿伦尼乌斯两点式包含ln(k2/k1)，因此对速率常数比值敏感。",
    },
]

DEMO_ANSWERS = [2, 2, 0, 1, 1, 1]


def run_game(interactive: bool) -> dict:
    details = []
    for index, item in enumerate(QUESTIONS, start=1):
        if interactive:
            print(f"\n{index}. {item['question']}")
            for i, option in enumerate(item["options"], start=1):
                print(f"  {i}. {option}")
            while True:
                raw = input("请输入选项编号：").strip()
                if raw.isdigit() and 1 <= int(raw) <= len(item["options"]):
                    answer = int(raw) - 1
                    break
                print("输入无效，请重新输入。")
        else:
            answer = DEMO_ANSWERS[index - 1]
        correct = answer == item["answer"]
        details.append(
            {
                "category": item["category"],
                "question": item["question"],
                "selected": item["options"][answer],
                "correct_answer": item["options"][item["answer"]],
                "correct": correct,
                "explanation": item["explanation"],
            }
        )
        if interactive:
            print("回答正确。" if correct else "回答不正确。")
            print(item["explanation"])

    score = sum(x["correct"] for x in details)
    category_scores = {}
    for category in sorted(set(x["category"] for x in details)):
        subset = [x for x in details if x["category"] == category]
        category_scores[category] = sum(x["correct"] for x in subset) / len(subset)
    return {
        "score": score,
        "total": len(QUESTIONS),
        "percentage": 100 * score / len(QUESTIONS),
        "category_scores": category_scores,
        "details": details,
    }


def save_outputs(result: dict) -> None:
    (OUTPUT_DIR / "知识游戏结果.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    lines = [
        f"总分：{result['score']}/{result['total']}（{result['percentage']:.1f}%）",
        "",
    ]
    for item in result["details"]:
        lines.extend(
            [
                item["question"],
                f"作答：{item['selected']}；正确答案：{item['correct_answer']}；{'正确' if item['correct'] else '错误'}",
                f"解析：{item['explanation']}",
                "",
            ]
        )
    (OUTPUT_DIR / "知识游戏结果说明.txt").write_text("\n".join(lines), encoding="utf-8")

    categories = list(result["category_scores"].keys())
    values = [100 * result["category_scores"][x] for x in categories]
    plt.figure(figsize=(8.8, 5.0))
    bars = plt.bar(categories, values)
    for bar, value in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value + 2,
            f"{value:.0f}%",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    plt.ylabel("正确率 / %")
    plt.title("知识游戏分类表现")
    plt.ylim(0, 110)
    plt.xticks(rotation=20, ha="right")
    plt.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "知识游戏分类表现.jpg", dpi=600, format="jpg", bbox_inches="tight"
    )
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="CAIA物理化学知识游戏")
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="启用命令行交互答题；默认运行演示答案",
    )
    args = parser.parse_args()
    result = run_game(args.interactive)
    save_outputs(result)
    print(
        f"知识游戏完成：{result['score']}/{result['total']}，输出目录：{OUTPUT_DIR.resolve()}"
    )


if __name__ == "__main__":
    main()
