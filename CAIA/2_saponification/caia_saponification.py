import os

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
plt.rcParams["mathtext.fontset"] = "stix"

OUTPUT_DIR = Path("outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

TIME_MIN = np.array([2, 4, 6, 8, 10, 15, 20, 25, 30, 35], dtype=float)
KAPPA_28 = np.array(
    [1895, 1727, 1592, 1490, 1408, 1264, 1153, 1076, 1020, 960], dtype=float
)
KAPPA_38 = np.array(
    [2100, 1891, 1690, 1576, 1522, 1375, 1295, 1240, 1197, 1167], dtype=float
)
CONCENTRATION_MOL_L = 0.0100
TEMPERATURES_C = [28.0, 38.0]
TEMPERATURES_K = [301.15, 311.15]

REPORT_FORMULAS = {
    28.0: {"B": 1522.42, "A": 0.0889386, "C": 601.712},
    38.0: {"B": 1575.49, "A": 0.183759, "C": 955.877},
}

USE_PYSR = False


def model(t: np.ndarray, B: float, A: float, C: float) -> np.ndarray:
    return B / (1.0 + A * np.asarray(t, dtype=float)) + C


def metrics(y_true: np.ndarray, y_pred: np.ndarray, t: np.ndarray) -> dict[str, float]:
    residual = np.asarray(y_true, float) - np.asarray(y_pred, float)
    sse = float(np.sum(residual**2))
    sst = float(np.sum((y_true - np.mean(y_true)) ** 2))
    idx = int(np.argmax(np.abs(residual)))
    return {
        "R2": 1.0 - sse / sst,
        "RMSE_uS_cm^-1": float(np.sqrt(np.mean(residual**2))),
        "MAE_uS_cm^-1": float(np.mean(np.abs(residual))),
        "SSE": sse,
        "max_abs_residual_uS_cm^-1": float(abs(residual[idx])),
        "max_abs_residual_time_min": float(t[idx]),
    }


def linear_fit(t: np.ndarray, kappa: np.ndarray, temperature_c: float) -> dict:
    f = REPORT_FORMULAS[temperature_c]
    kappa_0 = f["B"] + f["C"]
    kappa_inf = f["C"]
    z = (kappa_0 - kappa) / (kappa - kappa_inf)
    slope = float(np.dot(t, z) / np.dot(t, t))
    k = slope / CONCENTRATION_MOL_L
    prediction = model(t, kappa_0 - kappa_inf, CONCENTRATION_MOL_L * k, kappa_inf)
    residual = kappa - prediction
    transformed_sse = float(np.sum((z - slope * t) ** 2))
    transformed_sst = float(np.sum((z - np.mean(z)) ** 2))
    return {
        "method": "线性拟合",
        "temperature_C": temperature_c,
        "equation": "(κ0-κt)/(κt-κ∞)=c0·k·t",
        "parameters": {
            "slope_min^-1": slope,
            "k_L_mol^-1_min^-1": k,
            "kappa_0_uS_cm^-1": kappa_0,
            "kappa_inf_uS_cm^-1": kappa_inf,
            "linear_space_R2": 1.0 - transformed_sse / transformed_sst,
        },
        "prediction": prediction,
        "residual": residual,
        "metrics": metrics(kappa, prediction, t),
        "note": "κ0和κ∞采用实验报告中符号回归得到的端点参数，用于复现实验报告中的线性化处理。",
    }


def nonlinear_fit(t: np.ndarray, kappa: np.ndarray, temperature_c: float) -> dict:
    p0 = [float(kappa[0] - kappa[-1]), 0.1, float(kappa[-1] * 0.7)]
    popt, _ = curve_fit(
        model,
        t,
        kappa,
        p0=p0,
        bounds=([0.0, 1e-12, 0.0], [np.inf, np.inf, np.inf]),
        maxfev=100000,
    )
    B, A, C = map(float, popt)
    prediction = model(t, B, A, C)
    return {
        "method": "非线性拟合",
        "temperature_C": temperature_c,
        "equation": "κ(t)=B/(1+A·t)+C",
        "parameters": {
            "B_uS_cm^-1": B,
            "A_min^-1": A,
            "C_uS_cm^-1": C,
            "kappa_0_uS_cm^-1": B + C,
            "kappa_inf_uS_cm^-1": C,
            "k_L_mol^-1_min^-1": A / CONCENTRATION_MOL_L,
        },
        "prediction": prediction,
        "residual": kappa - prediction,
        "metrics": metrics(kappa, prediction, t),
        "note": "直接在原始电导率尺度上同时优化B、A和C。",
    }


def symbolic_fit_report_formula(
    t: np.ndarray, kappa: np.ndarray, temperature_c: float
) -> dict:
    f = REPORT_FORMULAS[temperature_c]
    prediction = model(t, f["B"], f["A"], f["C"])
    return {
        "method": "符号回归",
        "temperature_C": temperature_c,
        "equation": f"κ(t)={f['B']}/(1+{f['A']}·t)+{f['C']}",
        "parameters": {
            "B_uS_cm^-1": f["B"],
            "A_min^-1": f["A"],
            "C_uS_cm^-1": f["C"],
            "kappa_0_uS_cm^-1": f["B"] + f["C"],
            "kappa_inf_uS_cm^-1": f["C"],
            "k_L_mol^-1_min^-1": f["A"] / CONCENTRATION_MOL_L,
        },
        "prediction": prediction,
        "residual": kappa - prediction,
        "metrics": metrics(kappa, prediction, t),
        "note": "采用实验报告中已经获得的符号回归表达式，保证复算稳定。",
    }


def symbolic_fit_pysr(t: np.ndarray, kappa: np.ndarray, temperature_c: float) -> dict:
    try:
        from pysr import PySRRegressor
    except ImportError as exc:
        raise RuntimeError("未安装PySR，请先执行 pip install pysr") from exc
    reg = PySRRegressor(
        niterations=200,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["exp"],
        model_selection="best",
        parsimony=1e-3,
        random_state=42,
        deterministic=True,
        parallelism="serial",
        verbosity=0,
    )
    x = t.reshape(-1, 1)
    reg.fit(x, kappa)
    prediction = np.asarray(reg.predict(x), float)
    return {
        "method": "符号回归",
        "temperature_C": temperature_c,
        "equation": str(reg.sympy()),
        "parameters": {},
        "prediction": prediction,
        "residual": kappa - prediction,
        "metrics": metrics(kappa, prediction, t),
        "note": "由PySR自由搜索，应进一步检查量纲与物理意义。",
    }


def activation_energy(k1: float, k2: float, t1_k: float, t2_k: float) -> float:
    R = 8.314462618
    return R * math.log(k2 / k1) / (1.0 / t1_k - 1.0 / t2_k) / 1000.0


def to_jsonable(result: dict) -> dict:
    result = dict(result)
    result["prediction"] = np.asarray(result["prediction"]).tolist()
    result["residual"] = np.asarray(result["residual"]).tolist()
    return result


def run() -> tuple[dict[float, list[dict]], dict[str, float]]:
    datasets = {28.0: KAPPA_28, 38.0: KAPPA_38}
    all_results: dict[float, list[dict]] = {}
    for temp, values in datasets.items():
        methods = [
            linear_fit(TIME_MIN, values, temp),
            nonlinear_fit(TIME_MIN, values, temp),
        ]
        methods.append(
            symbolic_fit_pysr(TIME_MIN, values, temp)
            if USE_PYSR
            else symbolic_fit_report_formula(TIME_MIN, values, temp)
        )
        all_results[temp] = methods

    ea = {}
    for method in ["线性拟合", "非线性拟合", "符号回归"]:
        k28 = next(x for x in all_results[28.0] if x["method"] == method)[
            "parameters"
        ].get("k_L_mol^-1_min^-1")
        k38 = next(x for x in all_results[38.0] if x["method"] == method)[
            "parameters"
        ].get("k_L_mol^-1_min^-1")
        if k28 and k38:
            ea[method] = activation_energy(
                k28, k38, TEMPERATURES_K[0], TEMPERATURES_K[1]
            )
    return all_results, ea


def save_outputs(all_results: dict[float, list[dict]], ea: dict[str, float]) -> None:
    rows = []
    predictions = pd.DataFrame({"time_min": TIME_MIN})
    for temp, results in all_results.items():
        raw = KAPPA_28 if temp == 28.0 else KAPPA_38
        predictions[f"kappa_{int(temp)}C"] = raw
        for result in results:
            p = result["parameters"]
            rows.append(
                {
                    "温度_C": temp,
                    "方法": result["method"],
                    "方程": result["equation"],
                    "k_L_mol^-1_min^-1": p.get("k_L_mol^-1_min^-1", np.nan),
                    "活化能_kJ_mol^-1": ea.get(result["method"], np.nan),
                    **result["metrics"],
                    "说明": result["note"],
                }
            )
            predictions[f"{int(temp)}C_{result['method']}_prediction"] = result[
                "prediction"
            ]
            predictions[f"{int(temp)}C_{result['method']}_residual"] = result[
                "residual"
            ]

    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / "乙酸乙酯皂化_计算指标.csv", index=False, encoding="utf-8-sig"
    )
    predictions.to_csv(
        OUTPUT_DIR / "乙酸乙酯皂化_逐点预测.csv", index=False, encoding="utf-8-sig"
    )

    payload = {
        "project": "Chemistry Artificial Intelligence Assistant (CAIA)",
        "module": "乙酸乙酯皂化反应动力学分析",
        "conditions": {
            "concentration_mol_L": CONCENTRATION_MOL_L,
            "temperatures_C": TEMPERATURES_C,
        },
        "data": {
            "time_min": TIME_MIN.tolist(),
            "kappa_28C": KAPPA_28.tolist(),
            "kappa_38C": KAPPA_38.tolist(),
        },
        "activation_energy_kJ_mol": ea,
        "results": {
            str(temp): [to_jsonable(x) for x in results]
            for temp, results in all_results.items()
        },
    }
    (OUTPUT_DIR / "乙酸乙酯皂化_完整结果.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    lines = ["CAIA 乙酸乙酯皂化反应动力学分析", ""]
    for temp in TEMPERATURES_C:
        lines.append(f"【{temp:.1f} ℃】")
        for result in all_results[temp]:
            p = result["parameters"]
            m = result["metrics"]
            lines.extend(
                [
                    f"{result['method']}：k={p.get('k_L_mol^-1_min^-1', float('nan')):.8f} L·mol⁻¹·min⁻¹",
                    f"R²={m['R2']:.9f}，RMSE={m['RMSE_uS_cm^-1']:.6f} μS·cm⁻¹，最大绝对残差={m['max_abs_residual_uS_cm^-1']:.6f} μS·cm⁻¹",
                    f"说明：{result['note']}",
                ]
            )
        lines.append("")
    lines.append("【活化能】")
    for method, value in ea.items():
        lines.append(f"{method}：Ea={value:.6f} kJ·mol⁻¹")
    (OUTPUT_DIR / "乙酸乙酯皂化_结果说明.txt").write_text(
        "\n".join(lines), encoding="utf-8"
    )

    dense_t = np.linspace(TIME_MIN.min(), TIME_MIN.max(), 500)
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.2))
    for ax, temp, raw in zip(axes, TEMPERATURES_C, [KAPPA_28, KAPPA_38]):
        ax.scatter(TIME_MIN, raw, s=38, label="实验数据", zorder=5)
        for result in all_results[temp]:
            p = result["parameters"]
            if p:
                B = p.get(
                    "B_uS_cm^-1",
                    p.get("kappa_0_uS_cm^-1", 0) - p.get("kappa_inf_uS_cm^-1", 0),
                )
                A = p.get(
                    "A_min^-1", CONCENTRATION_MOL_L * p.get("k_L_mol^-1_min^-1", 0)
                )
                C = p.get("C_uS_cm^-1", p.get("kappa_inf_uS_cm^-1", 0))
                y_dense = model(dense_t, B, A, C)
            else:
                y_dense = np.interp(dense_t, TIME_MIN, result["prediction"])
            ax.plot(dense_t, y_dense, linewidth=1.7, label=result["method"])
        ax.set_title(f"{temp:.1f} ℃", fontsize=13)
        ax.set_xlabel("时间 t / min", fontsize=11)
        ax.set_ylabel(r"电导率 κ / (μS·cm$^{-1}$)", fontsize=11)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)
    fig.suptitle("乙酸乙酯皂化反应三种方法拟合结果", fontsize=14)
    fig.tight_layout()
    fig.savefig(
        OUTPUT_DIR / "乙酸乙酯皂化_三种方法拟合.jpg",
        dpi=600,
        format="jpg",
        bbox_inches="tight",
    )
    plt.close(fig)

    plt.figure(figsize=(7.4, 5.0))
    methods = list(ea.keys())
    values = [ea[m] for m in methods]
    bars = plt.bar(methods, values)
    for bar, value in zip(bars, values):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.25,
            f"{value:.2f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    plt.ylabel(r"活化能 Ea / (kJ·mol$^{-1}$)", fontsize=11)
    plt.title("不同数据处理方法得到的活化能", fontsize=14)
    plt.ylim(0, max(values) * 1.15)
    plt.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "乙酸乙酯皂化_活化能比较.jpg",
        dpi=600,
        format="jpg",
        bbox_inches="tight",
    )
    plt.close()


def main() -> None:
    all_results, ea = run()
    save_outputs(all_results, ea)
    print("乙酸乙酯皂化分析完成。输出目录：", OUTPUT_DIR.resolve())
    for temp, results in all_results.items():
        for result in results:
            k = result["parameters"].get("k_L_mol^-1_min^-1", float("nan"))
            print(
                f"{temp:.1f} ℃ {result['method']}: k={k:.8f}, R2={result['metrics']['R2']:.9f}"
            )
    for method, value in ea.items():
        print(f"{method}: Ea={value:.6f} kJ/mol")


if __name__ == "__main__":
    main()
