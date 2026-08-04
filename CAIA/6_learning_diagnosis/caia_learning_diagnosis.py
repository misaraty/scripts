import os

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

RECORDS = [
    ("S001", "反应级数判定", 6, 9, 10),
    ("S001", "积分速率方程", 5, 8, 10),
    ("S001", "理论实验关联", 5, 9, 10),
    ("S001", "数据处理", 6, 9, 10),
    ("S001", "误差分析", 4, 8, 10),
    ("S002", "反应级数判定", 7, 9, 10),
    ("S002", "积分速率方程", 6, 8, 10),
    ("S002", "理论实验关联", 5, 8, 10),
    ("S002", "数据处理", 5, 8, 10),
    ("S002", "误差分析", 4, 7, 10),
    ("S003", "反应级数判定", 5, 8, 10),
    ("S003", "积分速率方程", 4, 7, 10),
    ("S003", "理论实验关联", 4, 8, 10),
    ("S003", "数据处理", 5, 8, 10),
    ("S003", "误差分析", 3, 7, 10),
    ("S004", "反应级数判定", 8, 10, 10),
    ("S004", "积分速率方程", 7, 9, 10),
    ("S004", "理论实验关联", 6, 9, 10),
    ("S004", "数据处理", 7, 9, 10),
    ("S004", "误差分析", 5, 8, 10),
    ("S005", "反应级数判定", 6, 8, 10),
    ("S005", "积分速率方程", 5, 7, 10),
    ("S005", "理论实验关联", 4, 7, 10),
    ("S005", "数据处理", 4, 7, 10),
    ("S005", "误差分析", 3, 6, 10),
    ("S006", "反应级数判定", 7, 9, 10),
    ("S006", "积分速率方程", 5, 8, 10),
    ("S006", "理论实验关联", 5, 8, 10),
    ("S006", "数据处理", 6, 9, 10),
    ("S006", "误差分析", 4, 7, 10),
]


def analyze() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    df = pd.DataFrame(
        RECORDS,
        columns=[
            "student_id",
            "knowledge_point",
            "pre_score",
            "post_score",
            "max_score",
        ],
    )
    df["pre_mastery"] = df["pre_score"] / df["max_score"]
    df["post_mastery"] = df["post_score"] / df["max_score"]
    df["gain"] = df["post_score"] - df["pre_score"]
    denominator = df["max_score"] - df["pre_score"]
    df["normalized_gain"] = np.where(denominator > 0, df["gain"] / denominator, np.nan)

    knowledge = df.groupby("knowledge_point", as_index=False).agg(
        pre_mastery=("pre_mastery", "mean"),
        post_mastery=("post_mastery", "mean"),
        mean_gain=("gain", "mean"),
        mean_normalized_gain=("normalized_gain", "mean"),
    )
    knowledge["difficulty_rank"] = (
        knowledge["post_mastery"].rank(method="dense", ascending=True).astype(int)
    )

    students = df.groupby("student_id", as_index=False).agg(
        pre_mastery=("pre_mastery", "mean"),
        post_mastery=("post_mastery", "mean"),
        mean_gain=("gain", "mean"),
        mean_normalized_gain=("normalized_gain", "mean"),
    )
    students["recommendation"] = np.select(
        [students["post_mastery"] < 0.70, students["post_mastery"] < 0.85],
        ["基础巩固", "针对性提升"],
        default="综合迁移",
    )

    summary = {
        "record_count": int(len(df)),
        "student_count": int(df["student_id"].nunique()),
        "knowledge_point_count": int(df["knowledge_point"].nunique()),
        "pre_mean_score": float(df["pre_score"].mean()),
        "post_mean_score": float(df["post_score"].mean()),
        "mean_gain": float(df["gain"].mean()),
        "mean_normalized_gain": float(df["normalized_gain"].mean()),
        "weakest_knowledge_point": str(
            knowledge.sort_values("post_mastery").iloc[0]["knowledge_point"]
        ),
    }
    return df, knowledge, students, summary


def draw_heatmap(df: pd.DataFrame) -> None:
    matrix = df.pivot(
        index="student_id", columns="knowledge_point", values="post_mastery"
    )
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    image = ax.imshow(matrix.values, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(
        np.arange(len(matrix.columns)), labels=matrix.columns, rotation=25, ha="right"
    )
    ax.set_yticks(np.arange(len(matrix.index)), labels=matrix.index)
    ax.set_title("学生知识点掌握度热图", fontsize=14)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(
                j, i, f"{matrix.iloc[i, j]:.0%}", ha="center", va="center", fontsize=8
            )
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label("掌握度")
    fig.tight_layout()
    fig.savefig(
        OUTPUT_DIR / "学情掌握度热图.jpg", dpi=600, format="jpg", bbox_inches="tight"
    )
    plt.close(fig)


def draw_gain_chart(knowledge: pd.DataFrame) -> None:
    x = np.arange(len(knowledge))
    width = 0.35
    plt.figure(figsize=(9.4, 5.2))
    plt.bar(x - width / 2, knowledge["pre_mastery"] * 100, width, label="前测")
    plt.bar(x + width / 2, knowledge["post_mastery"] * 100, width, label="后测")
    plt.xticks(x, knowledge["knowledge_point"], rotation=20, ha="right")
    plt.ylabel("平均掌握度 / %", fontsize=11)
    plt.title("各知识点前后测掌握度", fontsize=14)
    plt.ylim(0, 105)
    plt.legend()
    plt.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "知识点前后测.jpg", dpi=600, format="jpg", bbox_inches="tight"
    )
    plt.close()


def main() -> None:
    df, knowledge, students, summary = analyze()
    df.to_csv(OUTPUT_DIR / "学情明细.csv", index=False, encoding="utf-8-sig")
    knowledge.to_csv(OUTPUT_DIR / "知识点诊断.csv", index=False, encoding="utf-8-sig")
    students.to_csv(
        OUTPUT_DIR / "学生诊断与建议.csv", index=False, encoding="utf-8-sig"
    )
    (OUTPUT_DIR / "学情汇总.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    text = (
        f"演示数据包含{summary['student_count']}名学生和{summary['knowledge_point_count']}个知识点，共{summary['record_count']}条前后测记录。\n"
        f"前测平均分为{summary['pre_mean_score']:.3f}/10，后测平均分为{summary['post_mean_score']:.3f}/10，平均增益为{summary['mean_gain']:.3f}。\n"
        f"按单条记录计算后再取平均的标准化增益为{summary['mean_normalized_gain']:.3f}。\n"
        f"当前平均掌握度最低的知识点是{summary['weakest_knowledge_point']}。"
    )
    (OUTPUT_DIR / "学情诊断说明.txt").write_text(text, encoding="utf-8")
    draw_heatmap(df)
    draw_gain_chart(knowledge)
    print(text)
    print("输出目录：", OUTPUT_DIR.resolve())


if __name__ == "__main__":
    main()
