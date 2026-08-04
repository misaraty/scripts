import os

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import math
import random
from pathlib import Path

import gradio as gr
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
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

# ==================== 动力学内置数据 ====================
SUC_T = np.array([5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55.0], float)
SUC_Y = np.array(
    [10.924, 9.0332, 7.414, 5.955, 4.741, 3.632, 2.675, 1.843, 1.098, 0.474, -0.083]
)
SAP_T = np.array([2, 4, 6, 8, 10, 15, 20, 25, 30, 35.0], float)
SAP_28 = np.array([1895, 1727, 1592, 1490, 1408, 1264, 1153, 1076, 1020, 960.0])
SAP_38 = np.array([2100, 1891, 1690, 1576, 1522, 1375, 1295, 1240, 1197, 1167.0])


def suc_model(t, amplitude, k, alpha_inf):
    return amplitude * np.exp(-k * t) + alpha_inf


def sap_model(t, amplitude, rate_factor, kappa_inf):
    return amplitude / (1 + rate_factor * t) + kappa_inf


def calculate_fit_metrics(observed, predicted):
    residual = observed - predicted
    sse = float(np.sum(residual**2))
    sst = float(np.sum((observed - observed.mean()) ** 2))
    r2 = 1 - sse / sst if sst > 0 else float("nan")
    rmse = float(np.sqrt(np.mean(residual**2)))
    max_residual = float(np.max(np.abs(residual)))
    return r2, rmse, max_residual


def analyze_kinetics(experiment):
    if experiment == "蔗糖水解 30.0 ℃":
        popt, _ = curve_fit(suc_model, SUC_T, SUC_Y, p0=[17, 0.027, -4], maxfev=100000)
        predicted = suc_model(SUC_T, *popt)
        r2, rmse, max_residual = calculate_fit_metrics(SUC_Y, predicted)
        _amplitude, k, alpha_inf = popt
        dense_t = np.linspace(SUC_T.min(), SUC_T.max(), 500)

        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        ax.scatter(SUC_T, SUC_Y, label="实验数据")
        ax.plot(dense_t, suc_model(dense_t, *popt), label="非线性拟合")
        ax.set_xlabel("时间 t / min")
        ax.set_ylabel("旋光度 α / °")
        ax.set_title("蔗糖水解非线性拟合")
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()

        text = (
            f"k = {k:.8f} min⁻¹\n"
            f"t₁/₂ = {math.log(2) / k:.6f} min\n"
            f"R² = {r2:.9f}\n"
            f"RMSE = {rmse:.6f}°\n"
            f"最大绝对残差 = {max_residual:.6f}°\n"
            f"α∞ = {alpha_inf:.6f}°"
        )
        return text, fig

    observed = SAP_28 if "28.0" in experiment else SAP_38
    temperature = 28.0 if "28.0" in experiment else 38.0
    popt, _ = curve_fit(
        sap_model,
        SAP_T,
        observed,
        p0=[1400, 0.1, 700],
        bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
        maxfev=100000,
    )
    predicted = sap_model(SAP_T, *popt)
    r2, rmse, max_residual = calculate_fit_metrics(observed, predicted)
    _amplitude, rate_factor, kappa_inf = popt
    k = rate_factor / 0.0100
    dense_t = np.linspace(SAP_T.min(), SAP_T.max(), 500)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.scatter(SAP_T, observed, label="实验数据")
    ax.plot(dense_t, sap_model(dense_t, *popt), label="非线性拟合")
    ax.set_xlabel("时间 t / min")
    ax.set_ylabel(r"电导率 κ / (μS·cm$^{-1}$)")
    ax.set_title(f"乙酸乙酯皂化 {temperature:.1f} ℃")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    text = (
        f"k = {k:.8f} L·mol⁻¹·min⁻¹\n"
        f"R² = {r2:.9f}\n"
        f"RMSE = {rmse:.6f} μS·cm⁻¹\n"
        f"最大绝对残差 = {max_residual:.6f} μS·cm⁻¹\n"
        f"κ∞ = {kappa_inf:.6f} μS·cm⁻¹"
    )
    return text, fig


# ==================== 简化的课程知识检索 ====================
DOCS = [
    (
        "一级反应",
        "一级反应满足-dc/dt=kc，积分后ln(c/c0)=-kt，半衰期t1/2=ln2/k，与初始浓度无关。",
    ),
    ("二级反应", "等初始浓度二级反应满足x/[a(a-x)]=kt，半衰期与初始浓度有关。"),
    (
        "蔗糖水解",
        "蔗糖水解可视为准一级反应，旋光度模型为alpha_t=alpha_inf+(alpha_0-alpha_inf)exp(-kt)。",
    ),
    (
        "皂化反应",
        "乙酸乙酯皂化可按二级反应处理，电导率模型为kappa_t=B/(1+At)+C，其中A=c0k。",
    ),
    (
        "模型评价",
        "不能只看R2，还应在原始数据尺度比较RMSE、MAE、最大绝对残差和参数物理意义。",
    ),
    (
        "可信人工智能",
        "确定性程序负责数值计算，课程知识库提供证据，大语言模型可用于解释，教师负责审核。",
    ),
]
VECTORIZER = TfidfVectorizer(analyzer="char", ngram_range=(2, 4))
DOC_MATRIX = VECTORIZER.fit_transform([f"{title} {text}" for title, text in DOCS])
REACTION_ORDERS = ("零级反应", "一级反应", "二级反应")


def detect_reaction_order(text):
    for order in REACTION_ORDERS:
        if order in text:
            return order
    return None


def rag_answer(question):
    question = question.strip()
    if not question:
        return "请输入问题。"

    query_vector = VECTORIZER.transform([question])
    base_scores = cosine_similarity(query_vector, DOC_MATRIX).ravel()
    adjusted_scores = base_scores.copy()
    target_order = detect_reaction_order(question)

    for index, (title, text) in enumerate(DOCS):
        document_text = f"{title} {text}"
        if title in question:
            adjusted_scores[index] += 0.20
        if target_order:
            if target_order in document_text:
                adjusted_scores[index] += 0.35
            elif any(
                order in document_text
                for order in REACTION_ORDERS
                if order != target_order
            ):
                adjusted_scores[index] -= 0.35

    adjusted_scores = np.clip(adjusted_scores, 0.0, 1.0)
    order = np.argsort(adjusted_scores)[::-1]
    result_limit = 1 if target_order else 2
    selected = [index for index in order if adjusted_scores[index] >= 0.05][
        :result_limit
    ]

    if not selected:
        return "当前演示知识库证据不足，请补充课程资料或由教师核查。"

    paragraphs = []
    for rank, index in enumerate(selected, start=1):
        title, text = DOCS[index]
        paragraphs.append(
            f"{text}（来源{rank}：{title}，检索得分{adjusted_scores[index]:.3f}）"
        )
    paragraphs.append("本页面只演示本地检索，不调用DeepSeek或其他大语言模型。")
    return "\n\n".join(paragraphs)


# ==================== 作业与报告反馈 ====================
def feedback(text, task_type):
    text = text.strip()
    if not text:
        return "请输入作业或报告内容。"

    suggestions = []
    if "R2" in text or "R²" in text:
        if not any(term in text for term in ["残差", "RMSE", "MAE", "物理意义"]):
            suggestions.append(
                "R²不能单独证明模型正确，还应比较残差、RMSE或参数的物理意义。"
            )

    if task_type == "实验报告":
        checks = [
            ("原始数据和单位", ["单位", "min", "μS"]),
            ("残差分析", ["残差", "RMSE"]),
            ("误差来源", ["误差", "恒温", "混合"]),
            ("结论中的活化能", ["结论", "活化能"]),
        ]
        missing = [
            label for label, words in checks if not any(word in text for word in words)
        ]
        if missing:
            suggestions.append("报告还可补充" + "、".join(missing) + "。")
    else:
        if "适用条件" not in text:
            suggestions.append("公式使用时应说明适用条件。")
        if "单位" not in text:
            suggestions.append("速率常数的单位需要核查并写明。")

    if not suggestions:
        suggestions.append("主要要素已经出现，可进一步压缩表述并说明参数的物理意义。")

    return (
        "本页面根据预设词语检查文本覆盖情况，只用于演示形成性反馈流程。\n\n"
        + "".join(suggestions)
        + "教师仍需结合原始数据和课程要求进行判断。"
    )


# ==================== 学情诊断 ====================
DEMO_MASTERY = {
    "反应级数判定": 0.86,
    "积分速率方程": 0.78,
    "理论实验关联": 0.82,
    "数据处理": 0.85,
    "误差分析": 0.72,
}


def diagnosis():
    weakest = min(DEMO_MASTERY, key=DEMO_MASTERY.get)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    labels = list(DEMO_MASTERY)
    values = [100 * DEMO_MASTERY[label] for label in labels]
    ax.bar(labels, values)
    ax.set_ylim(0, 100)
    ax.set_ylabel("掌握度 / %")
    ax.set_title("班级知识点掌握度（演示数据）")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    return (
        f"演示数据中平均掌握度最低的知识点是{weakest}，掌握度为{DEMO_MASTERY[weakest]:.0%}。"
        "可据此安排残差分析和实验误差练习。该结果不是实际教学效果证据。",
        fig,
    )


# ==================== 知识游戏 ====================
GAME = [
    ("一级反应半衰期与初始浓度？", ["成正比", "成反比", "无关"], "无关"),
    ("乙酸乙酯皂化主要用什么信号跟踪？", ["旋光度", "电导率", "压力"], "电导率"),
    ("曲线重合是否表示参数完全相同？", ["是", "否"], "否"),
]
CURRENT = {"item": GAME[0]}


def new_question():
    CURRENT["item"] = random.choice(GAME)
    question, options, _ = CURRENT["item"]
    return question, gr.update(choices=options, value=None), ""


def check_answer(answer):
    if not answer:
        return "请选择答案。"
    correct = CURRENT["item"][2]
    return "回答正确。" if answer == correct else f"回答不正确。正确答案是{correct}。"


def build_app():
    with gr.Blocks(title="Chemistry Artificial Intelligence Assistant (CAIA)") as demo:
        gr.Markdown("# Chemistry Artificial Intelligence Assistant (CAIA)")
        gr.Markdown("本页面为单文件功能演示，各模块采用简化实现。")

        with gr.Tab("动力学分析"):
            experiment = gr.Dropdown(
                ["蔗糖水解 30.0 ℃", "乙酸乙酯皂化 28.0 ℃", "乙酸乙酯皂化 38.0 ℃"],
                value="蔗糖水解 30.0 ℃",
                label="实验",
            )
            analyze_button = gr.Button("开始分析")
            result = gr.Textbox(label="计算结果", lines=8)
            plot = gr.Plot(label="拟合图")
            analyze_button.click(
                analyze_kinetics, inputs=experiment, outputs=[result, plot]
            )

        with gr.Tab("课程知识检索"):
            question = gr.Textbox(
                label="问题", value="一级反应半衰期与初始浓度有什么关系？"
            )
            rag_button = gr.Button("检索")
            rag_output = gr.Textbox(label="检索结果", lines=10)
            rag_button.click(rag_answer, inputs=question, outputs=rag_output)

        with gr.Tab("作业与报告反馈"):
            task_type = gr.Radio(
                ["理论作业", "实验报告"], value="实验报告", label="类型"
            )
            text = gr.Textbox(label="内容", lines=10)
            feedback_button = gr.Button("生成反馈")
            feedback_output = gr.Textbox(label="形成性反馈", lines=10)
            feedback_button.click(
                feedback, inputs=[text, task_type], outputs=feedback_output
            )

        with gr.Tab("学情诊断"):
            diagnosis_button = gr.Button("生成演示诊断")
            diagnosis_text = gr.Textbox(label="诊断")
            diagnosis_plot = gr.Plot(label="掌握度")
            diagnosis_button.click(diagnosis, outputs=[diagnosis_text, diagnosis_plot])

        with gr.Tab("知识游戏"):
            game_question = gr.Textbox(
                value=GAME[0][0], label="题目", interactive=False
            )
            choice = gr.Radio(GAME[0][1], label="选择")
            with gr.Row():
                new_button = gr.Button("换一题")
                check_button = gr.Button("提交")
            game_output = gr.Textbox(label="结果")
            new_button.click(new_question, outputs=[game_question, choice, game_output])
            check_button.click(check_answer, inputs=choice, outputs=game_output)

    return demo


def main():
    (OUTPUT_DIR / "启动说明.txt").write_text(
        "运行 caia_integrated_web_platform.py 后，浏览器访问 http://127.0.0.1:7860。设置CAIA_SHARE=1时会监听局域网地址，仅应在可信网络中使用。",
        encoding="utf-8",
    )

    lan_mode = os.getenv("CAIA_SHARE", "").strip() == "1"
    server_name = "0.0.0.0" if lan_mode else "127.0.0.1"

    username = os.getenv("CAIA_USERNAME", "").strip()
    password = os.getenv("CAIA_PASSWORD", "").strip()
    auth = (username, password) if username and password else None

    if lan_mode and auth is None:
        print(
            "提示：当前为局域网模式且未设置访问账号，请勿输入真实学生信息或其他敏感内容。"
        )

    build_app().launch(
        server_name=server_name,
        server_port=7860,
        auth=auth,
        share=False,
        show_error=False,
    )


if __name__ == "__main__":
    main()
