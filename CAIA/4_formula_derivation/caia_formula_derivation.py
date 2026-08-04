import os

os.chdir(os.path.split(os.path.realpath(__file__))[0])

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

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


def verify_formulas() -> dict:
    t, c0, k = sp.symbols("t c0 k", positive=True)
    zero = c0 - k * t
    first = c0 * sp.exp(-k * t)
    second = c0 / (1 + k * c0 * t)
    checks = {
        "零级反应": sp.simplify(sp.diff(zero, t) + k),
        "一级反应": sp.simplify(sp.diff(first, t) + k * first),
        "二级反应": sp.simplify(sp.diff(second, t) + k * second**2),
    }
    return {
        name: {"residual": str(value), "verified": bool(value == 0)}
        for name, value in checks.items()
    }


def build_markdown(checks: dict) -> str:
    return f"""# 简单级数反应公式推导与教学提示

## 1 零级反应
微分式：`-dc/dt=k`。积分得 `c=c0-kt`。作 `c-t` 图应为直线，斜率为 `-k`；半衰期为 `t1/2=c0/(2k)`。

## 2 一级反应
微分式：`-dc/dt=kc`。变量分离并积分得 `ln(c/c0)=-kt`，即 `c=c0 exp(-kt)`；半衰期为 `t1/2=ln2/k`，与初始浓度无关。

蔗糖水解在水和氢离子浓度近似不变时可视为准一级反应。旋光度模型为：

`alpha_t=alpha_inf+(alpha_0-alpha_inf)exp(-kt)`。

## 3 二级反应
等初始浓度条件下，微分式为 `dx/dt=k(a-x)^2`。积分得 `x/[a(a-x)]=kt`，半衰期为 `t1/2=1/(ka)`，与初始浓度有关。

乙酸乙酯皂化的电导率模型可写为：

`kappa_t=B/(1+At)+C`，其中 `A=c0 k`。

## 4 公式核验
- 零级反应：{checks['零级反应']['verified']}
- 一级反应：{checks['一级反应']['verified']}
- 二级反应：{checks['二级反应']['verified']}

## 5 分层练习
1. 基础：比较零级、一级和二级反应速率常数的单位。
2. 进阶：从一级反应积分式推导半衰期，并说明其为何与初始浓度无关。
3. 迁移：由乙酸乙酯皂化的二级反应积分式推导电导率与时间关系。
4. 评价：说明为什么曲线重合不代表模型参数完全一致。
"""


def draw_curves() -> None:
    t = np.linspace(0, 30, 500)
    c0 = 1.0
    zero = np.clip(c0 - 0.03 * t, 0, None)
    first = c0 * np.exp(-0.08 * t)
    second = c0 / (1 + 0.12 * c0 * t)
    plt.figure(figsize=(8.0, 5.2))
    plt.plot(t, zero, linewidth=2, label="零级反应")
    plt.plot(t, first, linewidth=2, label="一级反应")
    plt.plot(t, second, linewidth=2, label="二级反应")
    plt.xlabel("时间 t", fontsize=11)
    plt.ylabel(r"相对浓度 $c/c_0$", fontsize=11)
    plt.title("简单级数反应浓度随时间的变化", fontsize=14)
    plt.legend(fontsize=10)
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "简单级数反应曲线.jpg", dpi=600, format="jpg", bbox_inches="tight"
    )
    plt.close()


def main() -> None:
    checks = verify_formulas()
    (OUTPUT_DIR / "公式核验.json").write_text(
        json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (OUTPUT_DIR / "公式推导与练习.md").write_text(
        build_markdown(checks), encoding="utf-8"
    )
    draw_curves()
    print("公式推导与核验完成。输出目录：", OUTPUT_DIR.resolve())


if __name__ == "__main__":
    main()
