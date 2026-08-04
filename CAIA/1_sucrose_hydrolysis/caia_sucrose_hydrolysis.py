import os

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import linregress

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

TEMPERATURE_C = 30.0
HCL_CONCENTRATION_MOL_L = 2.0
SUCROSE_CONCENTRATION_MOL_L = 0.32
TIME_MIN = np.array([5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55], dtype=float)
OPTICAL_ROTATION_DEG = np.array(
    [10.924, 9.0332, 7.414, 5.955, 4.741, 3.632, 2.675, 1.843, 1.098, 0.474, -0.083],
    dtype=float,
)

REPORT_SYMBOLIC_A = 17.0230
REPORT_SYMBOLIC_K = 0.0270
REPORT_SYMBOLIC_C = -3.9557

USE_PYSR = False


def kinetic_model(
    t: np.ndarray, amplitude: float, k: float, alpha_inf: float
) -> np.ndarray:
    return amplitude * np.exp(-k * np.asarray(t, dtype=float)) + alpha_inf


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    residual = y_true - y_pred
    sse = float(np.sum(residual**2))
    sst = float(np.sum((y_true - np.mean(y_true)) ** 2))
    r2 = 1.0 - sse / sst if sst > 0 else float("nan")
    rmse = float(np.sqrt(np.mean(residual**2)))
    mae = float(np.mean(np.abs(residual)))
    max_index = int(np.argmax(np.abs(residual)))
    return {
        "R2": r2,
        "RMSE_deg": rmse,
        "MAE_deg": mae,
        "SSE": sse,
        "max_abs_residual_deg": float(abs(residual[max_index])),
        "max_abs_residual_time_min": float(TIME_MIN[max_index]),
    }


def linear_fit() -> dict:
    alpha_inf = REPORT_SYMBOLIC_C
    shifted = OPTICAL_ROTATION_DEG - alpha_inf
    if np.any(shifted <= 0):
        raise ValueError("线性拟合要求所有 α(t)-α∞ 大于0。")
    reg = linregress(TIME_MIN, np.log(shifted))
    k = -float(reg.slope)
    amplitude = float(np.exp(reg.intercept))
    prediction = kinetic_model(TIME_MIN, amplitude, k, alpha_inf)
    return {
        "method": "线性拟合",
        "equation": "ln[α(t)-α∞]=lnA-k·t",
        "parameters": {
            "A_deg": amplitude,
            "k_min^-1": k,
            "alpha_inf_deg": alpha_inf,
            "alpha_0_deg": amplitude + alpha_inf,
            "half_life_min": math.log(2) / k,
            "linear_space_R2": float(reg.rvalue**2),
        },
        "prediction": prediction,
        "residual": OPTICAL_ROTATION_DEG - prediction,
        "metrics": calculate_metrics(OPTICAL_ROTATION_DEG, prediction),
        "note": "α∞采用实验报告中符号回归式的常数项，用于复现实验报告中的线性化处理。",
    }


def nonlinear_fit() -> dict:
    p0 = [17.0, 0.027, -4.0]
    popt, _ = curve_fit(
        kinetic_model,
        TIME_MIN,
        OPTICAL_ROTATION_DEG,
        p0=p0,
        bounds=([-np.inf, 1e-12, -np.inf], [np.inf, np.inf, np.inf]),
        maxfev=100000,
    )
    amplitude, k, alpha_inf = map(float, popt)
    prediction = kinetic_model(TIME_MIN, amplitude, k, alpha_inf)
    return {
        "method": "非线性拟合",
        "equation": "α(t)=A·exp(-k·t)+C",
        "parameters": {
            "A_deg": amplitude,
            "k_min^-1": k,
            "alpha_inf_deg": alpha_inf,
            "alpha_0_deg": amplitude + alpha_inf,
            "half_life_min": math.log(2) / k,
        },
        "prediction": prediction,
        "residual": OPTICAL_ROTATION_DEG - prediction,
        "metrics": calculate_metrics(OPTICAL_ROTATION_DEG, prediction),
        "note": "直接在原始旋光度尺度上同时优化A、k和α∞。",
    }


def symbolic_regression_report_formula() -> dict:
    prediction = kinetic_model(
        TIME_MIN, REPORT_SYMBOLIC_A, REPORT_SYMBOLIC_K, REPORT_SYMBOLIC_C
    )
    return {
        "method": "符号回归",
        "equation": "α(t)=17.0230·exp(-0.0270·t)-3.9557",
        "parameters": {
            "A_deg": REPORT_SYMBOLIC_A,
            "k_min^-1": REPORT_SYMBOLIC_K,
            "alpha_inf_deg": REPORT_SYMBOLIC_C,
            "alpha_0_deg": REPORT_SYMBOLIC_A + REPORT_SYMBOLIC_C,
            "half_life_min": math.log(2) / REPORT_SYMBOLIC_K,
        },
        "prediction": prediction,
        "residual": OPTICAL_ROTATION_DEG - prediction,
        "metrics": calculate_metrics(OPTICAL_ROTATION_DEG, prediction),
        "note": "默认使用实验报告中已经获得的符号回归显式表达式，保证复算结果稳定。",
    }


def symbolic_regression_pysr() -> dict:
    try:
        from pysr import PySRRegressor
    except ImportError as exc:
        raise RuntimeError("未安装PySR，请先执行 pip install pysr") from exc

    model = PySRRegressor(
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
    x = TIME_MIN.reshape(-1, 1)
    model.fit(x, OPTICAL_ROTATION_DEG)
    prediction = np.asarray(model.predict(x), dtype=float)
    equation = str(model.sympy())
    return {
        "method": "符号回归",
        "equation": equation,
        "parameters": {},
        "prediction": prediction,
        "residual": OPTICAL_ROTATION_DEG - prediction,
        "metrics": calculate_metrics(OPTICAL_ROTATION_DEG, prediction),
        "note": "由PySR自由搜索得到，应进一步检查量纲、边界和物理意义。",
    }


def serializable(result: dict) -> dict:
    output = dict(result)
    output["prediction"] = np.asarray(result["prediction"]).tolist()
    output["residual"] = np.asarray(result["residual"]).tolist()
    return output


def save_outputs(results: list[dict]) -> None:
    prediction_df = pd.DataFrame(
        {
            "time_min": TIME_MIN,
            "optical_rotation_deg": OPTICAL_ROTATION_DEG,
        }
    )
    for result in results:
        prediction_df[f"{result['method']}_prediction"] = result["prediction"]
        prediction_df[f"{result['method']}_residual"] = result["residual"]
    prediction_df.to_csv(
        OUTPUT_DIR / "蔗糖水解_逐点预测.csv", index=False, encoding="utf-8-sig"
    )

    metric_rows = []
    for result in results:
        row = {
            "方法": result["method"],
            "方程": result["equation"],
            "k_min^-1": result["parameters"].get("k_min^-1", np.nan),
            "半衰期_min": result["parameters"].get("half_life_min", np.nan),
            **result["metrics"],
            "说明": result["note"],
        }
        metric_rows.append(row)
    pd.DataFrame(metric_rows).to_csv(
        OUTPUT_DIR / "蔗糖水解_计算指标.csv", index=False, encoding="utf-8-sig"
    )

    payload = {
        "project": "Chemistry Artificial Intelligence Assistant (CAIA)",
        "module": "蔗糖水解反应动力学分析",
        "conditions": {
            "temperature_C": TEMPERATURE_C,
            "HCl_concentration_mol_L": HCL_CONCENTRATION_MOL_L,
            "sucrose_concentration_mol_L": SUCROSE_CONCENTRATION_MOL_L,
        },
        "data": {
            "time_min": TIME_MIN.tolist(),
            "optical_rotation_deg": OPTICAL_ROTATION_DEG.tolist(),
        },
        "results": [serializable(result) for result in results],
    }
    (OUTPUT_DIR / "蔗糖水解_完整结果.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    lines = [
        "CAIA 蔗糖水解反应动力学分析",
        f"实验条件：{TEMPERATURE_C:.1f} ℃，盐酸浓度 {HCL_CONCENTRATION_MOL_L:.1f} mol/L，蔗糖初始浓度 {SUCROSE_CONCENTRATION_MOL_L:.2f} mol/L。",
        "",
    ]
    for result in results:
        p = result["parameters"]
        m = result["metrics"]
        lines.extend(
            [
                f"【{result['method']}】",
                f"方程：{result['equation']}",
                f"k = {p.get('k_min^-1', float('nan')):.8f} min⁻¹",
                f"t1/2 = {p.get('half_life_min', float('nan')):.6f} min",
                f"R² = {m['R2']:.9f}",
                f"RMSE = {m['RMSE_deg']:.6f}°，MAE = {m['MAE_deg']:.6f}°",
                f"最大绝对残差 = {m['max_abs_residual_deg']:.6f}°，出现在 t = {m['max_abs_residual_time_min']:.0f} min",
                f"说明：{result['note']}",
                "",
            ]
        )
    (OUTPUT_DIR / "蔗糖水解_结果说明.txt").write_text(
        "\n".join(lines), encoding="utf-8"
    )

    dense_t = np.linspace(TIME_MIN.min(), TIME_MIN.max(), 500)
    plt.figure(figsize=(8.2, 5.4))
    plt.scatter(
        TIME_MIN, OPTICAL_ROTATION_DEG, s=42, marker="o", label="实验数据", zorder=5
    )
    for result in results:
        if result["method"] == "线性拟合":
            p = result["parameters"]
            dense_y = kinetic_model(
                dense_t, p["A_deg"], p["k_min^-1"], p["alpha_inf_deg"]
            )
        elif result["method"] == "非线性拟合":
            p = result["parameters"]
            dense_y = kinetic_model(
                dense_t, p["A_deg"], p["k_min^-1"], p["alpha_inf_deg"]
            )
        else:
            if result["parameters"]:
                p = result["parameters"]
                dense_y = kinetic_model(
                    dense_t, p["A_deg"], p["k_min^-1"], p["alpha_inf_deg"]
                )
            else:
                dense_y = np.interp(dense_t, TIME_MIN, result["prediction"])
        plt.plot(dense_t, dense_y, linewidth=1.8, label=result["method"])
    plt.xlabel("时间 t / min", fontsize=12)
    plt.ylabel("旋光度 α / °", fontsize=12)
    plt.title("蔗糖水解反应三种方法拟合结果", fontsize=14)
    plt.legend(frameon=True, fontsize=10)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "蔗糖水解_三种方法拟合.jpg",
        dpi=600,
        format="jpg",
        bbox_inches="tight",
    )
    plt.close()

    plt.figure(figsize=(8.2, 5.4))
    for result in results:
        plt.plot(
            TIME_MIN,
            result["residual"],
            marker="o",
            linewidth=1.4,
            label=result["method"],
        )
    plt.axhline(0.0, linewidth=1.0)
    plt.xlabel("时间 t / min", fontsize=12)
    plt.ylabel("残差 / °", fontsize=12)
    plt.title("蔗糖水解反应残差比较", fontsize=14)
    plt.legend(frameon=True, fontsize=10)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "蔗糖水解_残差比较.jpg", dpi=600, format="jpg", bbox_inches="tight"
    )
    plt.close()


def main() -> None:
    results = [linear_fit(), nonlinear_fit()]
    results.append(
        symbolic_regression_pysr() if USE_PYSR else symbolic_regression_report_formula()
    )
    save_outputs(results)
    print("蔗糖水解分析完成。输出目录：", OUTPUT_DIR.resolve())
    for result in results:
        k = result["parameters"].get("k_min^-1", float("nan"))
        print(f"{result['method']}: k={k:.8f} min^-1, R2={result['metrics']['R2']:.9f}")


if __name__ == "__main__":
    main()
