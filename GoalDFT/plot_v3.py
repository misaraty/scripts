import os, csv, math
from pathlib import Path

os.chdir(os.path.split(os.path.realpath(__file__))[0])
import numpy as np
import matplotlib as mpl

mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

TAB = mpl.colormaps["tab10"].colors
COLORS = {
    "blue": TAB[0],
    "orange": TAB[1],
    "green": TAB[2],
    "red": TAB[3],
    "purple": TAB[4],
    "brown": TAB[5],
    "pink": TAB[6],
    "gray": TAB[7],
    "olive": TAB[8],
    "cyan": TAB[9],
    "lightgray": (0.88, 0.88, 0.88),
    "black": (0.08, 0.08, 0.08),
}
COLOR_ORDER = [TAB[i] for i in range(10)]
FONT_FAMILY = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
AXIS_LABEL_SIZE = 13.5
TICK_LABEL_SIZE = 11.5
LEGEND_SIZE = 10.5
PANEL_TITLE_SIZE = 12.0
PANEL_LABEL_SIZE = 12.5
ANNOTATION_SIZE = 10.0
TABLE_TITLE_SIZE = 11.5
TABLE_FONT_SIZE = 9.0
LINE_WIDTH = 2.0
MARKER_SIZE = 5.5
SPINE_WIDTH = 0.9
PANEL_LABEL_X = -0.155
PANEL_LABEL_Y = 1.015
DPI = 600
SAFETY_FACTOR = 0.65
NEAR_THRESHOLD_FACTOR = 1.10

mpl.rcParams["font.sans-serif"] = FONT_FAMILY
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["mathtext.fontset"] = "custom"
mpl.rcParams["mathtext.rm"] = "Arial"
mpl.rcParams["mathtext.it"] = "Arial:italic"
mpl.rcParams["mathtext.bf"] = "Arial:bold"
mpl.rcParams["mathtext.default"] = "it"
mpl.rcParams["axes.unicode_minus"] = False
mpl.rcParams["axes.labelsize"] = AXIS_LABEL_SIZE
mpl.rcParams["xtick.labelsize"] = TICK_LABEL_SIZE
mpl.rcParams["ytick.labelsize"] = TICK_LABEL_SIZE
mpl.rcParams["legend.fontsize"] = LEGEND_SIZE
mpl.rcParams["axes.edgecolor"] = COLORS["black"]
mpl.rcParams["axes.labelcolor"] = COLORS["black"]
mpl.rcParams["xtick.color"] = COLORS["black"]
mpl.rcParams["ytick.color"] = COLORS["black"]
mpl.rcParams["text.color"] = COLORS["black"]
mpl.rcParams["axes.linewidth"] = SPINE_WIDTH
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"
mpl.rcParams["xtick.major.width"] = 0.8
mpl.rcParams["ytick.major.width"] = 0.8
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"

HERE = Path(os.path.dirname(os.path.realpath(__file__)))
DATA_ROOT = HERE / "data"
if not DATA_ROOT.exists():
    DATA_ROOT = HERE.parent / "data"
if not DATA_ROOT.exists():
    raise FileNotFoundError(
        "Cannot find data/main and data/si beside plot_v1.py or one directory above it."
    )
MAIN = str(DATA_ROOT / "main")
SI = str(DATA_ROOT / "si")
FIG_MAIN = os.path.join("figures", "main")
FIG_SI = os.path.join("figures", "si")
TAB_MAIN = os.path.join("tables", "main")
TAB_SI = os.path.join("tables", "si")
for d in [FIG_MAIN, FIG_SI, TAB_MAIN, TAB_SI]:
    os.makedirs(d, exist_ok=True)
    for name in os.listdir(d):
        if name.lower().endswith((".jpg", ".jpeg", ".png", ".tsv")):
            os.remove(os.path.join(d, name))

SYSTEM_LABELS = {
    "Si": "Si",
    "C": "C",
    "MgO": "MgO",
    "LiF": "LiF",
    "Si_displaced": "Si",
    "MgO_displaced": "MgO",
    "TiO2_displaced": r"TiO$_2$",
}
SYSTEM_ORDER_ENERGY = ["Si", "C", "MgO", "LiF"]
SYSTEM_ORDER_FORCE = ["Si_displaced", "MgO_displaced", "TiO2_displaced"]


def read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fval(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def ival(x):
    try:
        return int(float(x))
    except Exception:
        return None


def bval(x):
    return str(x).lower() in ("true", "1", "yes")


def estimate_val(r):
    return fval(r.get("estimated_controllable_error", r.get("total_estimate", "nan")))


def syslabel(name):
    return SYSTEM_LABELS.get(name, name.replace("_displaced", ""))


def savefig(fig, path):
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def fixed(x, digits=2):
    if not math.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def sci(x, digits=2):
    if not math.isfinite(x):
        return "NA"
    if x == 0:
        return "0"
    return f"{x:.{digits}e}"


def time_text(x):
    if not math.isfinite(x):
        return "NA"
    return f"{x:.2f}" if x < 100 else f"{x:.1f}"


def geometric_mean(vals):
    vals = [v for v in vals if v > 0 and math.isfinite(v)]
    return (
        math.exp(sum(math.log(v) for v in vals) / len(vals)) if vals else float("nan")
    )


def style_axis(ax, grid="both"):
    for side in ["left", "bottom", "top", "right"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(SPINE_WIDTH)
        ax.spines[side].set_edgecolor(COLORS["black"])
    ax.tick_params(length=3.5, pad=3, top=False, right=False)
    if grid == "both":
        ax.grid(
            True, which="major", color=COLORS["lightgray"], linewidth=0.75, zorder=0
        )
    elif grid == "y":
        ax.grid(
            True,
            axis="y",
            which="major",
            color=COLORS["lightgray"],
            linewidth=0.75,
            zorder=0,
        )
    ax.set_axisbelow(True)


def style_twin_axis(ax, color=COLORS["black"]):
    ax.spines["top"].set_visible(True)
    ax.spines["top"].set_linewidth(SPINE_WIDTH)
    ax.spines["top"].set_edgecolor(COLORS["black"])
    ax.spines["right"].set_visible(True)
    ax.spines["right"].set_linewidth(SPINE_WIDTH)
    ax.spines["right"].set_edgecolor(COLORS["black"])
    ax.tick_params(
        axis="y", colors=color, labelsize=TICK_LABEL_SIZE, right=True, length=3.5, pad=3
    )
    ax.tick_params(axis="x", top=False)


def panel_label(ax, letter):
    ax.text(
        PANEL_LABEL_X,
        PANEL_LABEL_Y,
        f"({letter})",
        transform=ax.transAxes,
        fontsize=PANEL_LABEL_SIZE,
        fontweight="normal",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def write_table(tsv_columns, display_columns, rows, tsv_path, jpg_path, title):
    with open(tsv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(tsv_columns)
        w.writerows(rows)
    weights = []
    for c, header in enumerate(display_columns):
        header_len = max(len(part) for part in str(header).replace("$", "").split("\n"))
        body_len = max([len(str(r[c])) for r in rows] + [0])
        weights.append(max(7.0, min(25.0, max(header_len, body_len) + 2.0)))
    weights = np.asarray(weights, dtype=float)
    col_widths = 0.98 * weights / weights.sum()
    fig_w = max(9.0, min(16.5, 0.12 * float(weights.sum())))
    fig_h = max(2.2, 0.43 * (len(rows) + 1) + 1.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    ax.set_title(title, fontsize=TABLE_TITLE_SIZE, loc="left", pad=8)
    t = ax.table(
        cellText=rows,
        colLabels=display_columns,
        colWidths=col_widths,
        loc="center",
        cellLoc="center",
        colLoc="center",
    )
    t.auto_set_font_size(False)
    t.set_fontsize(TABLE_FONT_SIZE)
    t.scale(1.0, 1.24)
    nrows = len(rows)
    for (r, c), cell in t.get_celld().items():
        cell.set_facecolor("white")
        cell.set_edgecolor(COLORS["black"])
        if r == 0:
            cell.visible_edges = "TB"
            cell.set_linewidth(0.9)
            cell.set_text_props(weight="bold")
            cell.set_height(cell.get_height() * 1.45)
        elif r == nrows:
            cell.visible_edges = "B"
            cell.set_linewidth(0.8)
        else:
            cell.visible_edges = ""
        if c == 0 and r > 0:
            cell.get_text().set_ha("left")
    savefig(fig, jpg_path)


def fig2_reliability():
    e = read_tsv(os.path.join(MAIN, "fig5_target_reliability.dat"))
    f = read_tsv(os.path.join(SI, "figS3_force_cost.dat"))
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.55))
    for j, (ax, data, order, target_key, error_key, title, xlabel, ylabel) in enumerate(
        [
            (
                axes[0],
                e,
                SYSTEM_ORDER_ENERGY,
                "target_meV_per_atom",
                "actual_error_meV_per_atom",
                "Energy",
                r"Requested tolerance, $\tau_E$ (meV $\cdot$ atom$^{-1}$)",
                r"Realized error, $e_E$ (meV $\cdot$ atom$^{-1}$)",
            ),
            (
                axes[1],
                f,
                SYSTEM_ORDER_FORCE,
                "target_eV_per_A",
                "actual_error_eV_per_A",
                "Force",
                r"Requested tolerance, $\tau_F$ (eV $\cdot$ Å$^{-1}$)",
                r"Realized error, $e_F$ (eV $\cdot$ Å$^{-1}$)",
            ),
        ]
    ):
        vals = []
        for i, system in enumerate(order):
            rr = [r for r in data if r["system"] == system]
            x = [fval(r[target_key]) for r in rr]
            y = [fval(r[error_key]) for r in rr]
            vals += x + y
            ax.scatter(
                x, y, s=44, color=COLOR_ORDER[i], label=syslabel(system), zorder=3
            )
        lim = max(vals) * 1.08
        ax.plot(
            [0, lim],
            [0, lim],
            "--",
            color=COLORS["black"],
            linewidth=1.0,
            label=r"$e_Q=\tau_Q$",
        )
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=PANEL_TITLE_SIZE)
        ax.legend(frameon=False, handletextpad=0.4, labelspacing=0.35)
        style_axis(ax)
        panel_label(ax, chr(ord("a") + j))
    fig.subplots_adjust(wspace=0.30)
    savefig(fig, os.path.join(FIG_MAIN, "Fig2_target_reliability.jpg"))


def fig3_cost():
    e = read_tsv(os.path.join(MAIN, "fig2_cost_to_target.dat"))
    f = read_tsv(os.path.join(SI, "figS3_force_cost.dat"))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.55))
    energy_targets = [10.0, 5.0, 2.0]
    x = np.arange(len(SYSTEM_ORDER_ENERGY))
    w = 0.24
    for j, target in enumerate(energy_targets):
        vals = [
            next(
                fval(r["uniform_over_adaptive_speedup"])
                for r in e
                if r["system"] == system
                and abs(fval(r["target_meV_per_atom"]) - target) < 1e-12
            )
            for system in SYSTEM_ORDER_ENERGY
        ]
        axes[0].bar(
            x + (j - 1) * w,
            vals,
            width=w,
            color=COLOR_ORDER[j],
            label=rf"$\tau_E={target:g}$ meV $\cdot$ atom$^{{-1}}$",
        )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([syslabel(s) for s in SYSTEM_ORDER_ENERGY])
    axes[0].set_title("Energy targets", fontsize=PANEL_TITLE_SIZE)
    axes[0].legend(frameon=False)

    force_targets = [0.05, 0.02]
    x2 = np.arange(len(SYSTEM_ORDER_FORCE))
    w2 = 0.32
    for j, target in enumerate(force_targets):
        vals = [
            next(
                fval(r["uniform_over_adaptive_speedup"])
                for r in f
                if r["system"] == system
                and abs(fval(r["target_eV_per_A"]) - target) < 1e-12
            )
            for system in SYSTEM_ORDER_FORCE
        ]
        axes[1].bar(
            x2 + (j - 0.5) * w2,
            vals,
            width=w2,
            color=COLOR_ORDER[j],
            label=rf"$\tau_F={target:g}$ eV $\cdot$ Å$^{{-1}}$",
        )
    axes[1].set_xticks(x2)
    axes[1].set_xticklabels([syslabel(s) for s in SYSTEM_ORDER_FORCE])
    axes[1].set_title("Force targets", fontsize=PANEL_TITLE_SIZE)
    axes[1].legend(frameon=False)

    for j, ax in enumerate(axes):
        ax.axhline(1.0, linestyle="--", color=COLORS["black"], linewidth=1.0)
        ax.set_ylabel(r"Wall-time ratio, $t_{\mathrm{uniform}}/t_{\mathrm{GoalDFT}}$")
        ax.text(
            0.98,
            0.56,
            "GoalDFT faster",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=ANNOTATION_SIZE,
            color=COLORS["gray"],
        )
        ax.text(
            0.98,
            0.48,
            "Uniform faster",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=ANNOTATION_SIZE,
            color=COLORS["gray"],
        )
        style_axis(ax, grid="y")
        panel_label(ax, chr(ord("a") + j))
    fig.subplots_adjust(wspace=0.30)
    savefig(fig, os.path.join(FIG_MAIN, "Fig3_efficiency_vs_uniform.jpg"))


def fig4_calibration():
    e = read_tsv(os.path.join(MAIN, "fig1_estimator_calibration.dat"))
    f = read_tsv(os.path.join(SI, "figS2_force_estimator_calibration.dat"))
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.75))
    specs = [
        (axes[0], e, SYSTEM_ORDER_ENERGY, "Energy", r"meV $\cdot$ atom$^{-1}$", "E"),
        (axes[1], f, SYSTEM_ORDER_FORCE, "Force", r"eV $\cdot$ Å$^{-1}$", "F"),
    ]
    for j, (ax, rows, order, title, unit, qsym) in enumerate(specs):
        allv = []
        for i, system in enumerate(order):
            rr = [
                r
                for r in rows
                if r["system"] == system
                and estimate_val(r) > 0
                and fval(r["actual_error"]) > 0
            ]
            x = [estimate_val(r) for r in rr]
            y = [fval(r["actual_error"]) for r in rr]
            allv += x + y
            ax.scatter(
                x, y, s=32, color=COLOR_ORDER[i], label=syslabel(system), zorder=3
            )
        lo = min(allv)
        hi = max(allv)
        grid = np.logspace(np.log10(lo), np.log10(hi), 300)
        ax.fill_between(
            grid, 0.5 * grid, 2.0 * grid, color=COLORS["lightgray"], zorder=0
        )
        ax.plot(
            [lo, hi],
            [lo, hi],
            "--",
            color=COLORS["black"],
            linewidth=1.0,
            label=r"$e_Q=\widehat e_Q$",
        )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(rf"Estimated controllable error, $\widehat e_{qsym}$ ({unit})")
        ax.set_ylabel(rf"Realized total error, $e_{qsym}$ ({unit})")
        ax.set_title(title, fontsize=PANEL_TITLE_SIZE)
        ax.legend(frameon=False, handletextpad=0.4, labelspacing=0.35)
        ax.text(
            0.04,
            0.94,
            "gray band: factor of 2",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=ANNOTATION_SIZE,
            color=COLORS["gray"],
        )
        style_axis(ax)
        panel_label(ax, chr(ord("a") + j))
    fig.subplots_adjust(wspace=0.36)
    savefig(fig, os.path.join(FIG_MAIN, "Fig4_estimator_calibration.jpg"))


def fig5_budget_paths():
    e = read_tsv(os.path.join(MAIN, "fig3_error_budget.dat"))
    f = read_tsv(os.path.join(SI, "figS2_force_estimator_calibration.dat"))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.55))
    specs = [
        (
            e,
            "Si",
            "energy",
            5.0,
            r"Si: $\tau_E=5$ meV $\cdot$ atom$^{-1}$",
            r"Error (meV $\cdot$ atom$^{-1}$)",
        ),
        (
            f,
            "TiO2_displaced",
            "force",
            0.02,
            r"TiO$_2$: $\tau_F=0.02$ eV $\cdot$ Å$^{-1}$",
            r"Error (eV $\cdot$ Å$^{-1}$)",
        ),
    ]
    for j, (ax, (src, system, goal, target, title, ylabel)) in enumerate(
        zip(axes, specs)
    ):
        rr = [
            r
            for r in src
            if r["system"] == system
            and r["goal"] == goal
            and abs(fval(r["target"]) - target) < 1e-12
        ]
        rr.sort(key=lambda r: ival(r["step"]))
        steps = [ival(r["step"]) for r in rr]
        series = [
            ("basis_estimate", COLORS["blue"], r"$\eta_{\mathrm{PW}}$"),
            ("kpoint_estimate", COLORS["orange"], r"$\eta_k$"),
            ("scf_estimate", COLORS["green"], r"$\eta_{\mathrm{SCF}}$"),
        ]
        for key, c, label in series:
            ax.plot(
                steps,
                [fval(r[key]) for r in rr],
                marker="o",
                markersize=MARKER_SIZE,
                linewidth=LINE_WIDTH,
                color=c,
                label=label,
            )
        ax.plot(
            steps,
            [estimate_val(r) for r in rr],
            marker="o",
            markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH,
            color=COLORS["purple"],
            label=r"$\widehat e_Q$",
        )
        ax.plot(
            steps,
            [fval(r["actual_error"]) for r in rr],
            marker="s",
            markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH,
            color=COLORS["red"],
            label=r"$e_Q$",
        )
        ax.axhline(
            SAFETY_FACTOR * target,
            linestyle=":",
            color=COLORS["gray"],
            linewidth=1.0,
            label=r"$s\tau_Q$",
        )
        ax.set_yscale("log")
        ax.set_xlabel("Control step")
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=PANEL_TITLE_SIZE)
        style_axis(ax)
        panel_label(ax, chr(ord("a") + j))
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(
        h,
        l,
        ncol=6,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        columnspacing=1.0,
        handletextpad=0.4,
    )
    fig.subplots_adjust(wspace=0.30, top=0.80)
    savefig(fig, os.path.join(FIG_MAIN, "Fig5_error_budget_paths.jpg"))


def table_main():
    e = read_tsv(os.path.join(MAIN, "fig2_cost_to_target.dat"))
    f = read_tsv(os.path.join(SI, "figS3_force_cost.dat"))
    rows = []
    for goal, data, target_key, error_key in [
        ("Energy", e, "target_meV_per_atom", "actual_error_meV_per_atom"),
        ("Force", f, "target_eV_per_A", "actual_error_eV_per_A"),
    ]:
        ratios = [fval(r[error_key]) / fval(r[target_key]) for r in data]
        speeds = [fval(r["uniform_over_adaptive_speedup"]) for r in data]
        rows.append(
            [
                goal,
                str(len(data)),
                f"{sum(bval(r['success']) for r in data)}/{len(data)}",
                fixed(float(np.median(ratios)), 2),
                fixed(max(ratios), 2),
                fixed(geometric_mean(speeds), 2),
                fixed(max(speeds), 2),
                fixed(min(speeds), 2),
            ]
        )
    tsv_cols = [
        "Goal",
        "Cases",
        "Target success",
        "Median error/target",
        "Worst error/target",
        "Geometric mean uniform/GoalDFT wall-time ratio",
        "Best uniform/GoalDFT ratio",
        "Worst uniform/GoalDFT ratio",
    ]
    display_cols = [
        "Goal",
        "Cases",
        "Target\nsuccess",
        "Median\n" + r"$e_Q/\tau_Q$",
        "Worst\n" + r"$e_Q/\tau_Q$",
        "Geometric mean wall-time ratio\n"
        + r"$t_{\mathrm{uniform}}/t_{\mathrm{GoalDFT}}$",
        "Best\nratio",
        "Worst\nratio",
    ]
    write_table(
        tsv_cols,
        display_cols,
        rows,
        os.path.join(TAB_MAIN, "Table1_summary.tsv"),
        os.path.join(TAB_MAIN, "Table1_summary.jpg"),
        "Table 1. Reliability and computational efficiency of GoalDFT.",
    )


def figs_si():
    rows = read_tsv(os.path.join(MAIN, "fig4_parameter_path.dat"))
    specs = [("Si", 10.0), ("C", 5.0), ("MgO", 2.0)]
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.15))
    for j, (ax, (system, target)) in enumerate(zip(axes, specs)):
        rr = [
            r
            for r in rows
            if r["system"] == system and abs(fval(r["target"]) - target) < 1e-12
        ]
        rr.sort(key=lambda r: ival(r["step"]))
        steps = [ival(r["step"]) for r in rr]
        (line1,) = ax.plot(
            steps,
            [fval(r["Ecut_Ha"]) for r in rr],
            marker="o",
            linewidth=LINE_WIDTH,
            color=COLORS["blue"],
            label=r"$E_{\mathrm{cut}}$",
        )
        ax2 = ax.twinx()
        (line2,) = ax2.plot(
            steps,
            [fval(r["kgrid_n"]) for r in rr],
            marker="s",
            linewidth=LINE_WIDTH,
            color=COLORS["orange"],
            label=r"$n_k$",
        )
        ax.set_xlabel("Control step")
        ax.set_ylabel(
            r"Plane-wave cutoff, $E_{\mathrm{cut}}$ (Ha)", color=COLORS["blue"]
        )
        ax2.set_ylabel(r"$k$-grid size, $n_k$", color=COLORS["orange"])
        style_twin_axis(ax2, COLORS["orange"])
        ax.set_title(
            rf"{syslabel(system)}, $\tau_E={target:g}$ meV $\cdot$ atom$^{{-1}}$",
            fontsize=PANEL_TITLE_SIZE,
        )
        style_axis(ax)
        panel_label(ax, chr(ord("a") + j))
        if j == 0:
            ax.legend(
                [line1, line2],
                [r"$E_{\mathrm{cut}}$", r"$n_k$"],
                frameon=False,
                loc="best",
            )
    fig.subplots_adjust(wspace=0.58)
    savefig(fig, os.path.join(FIG_SI, "FigS1_parameter_paths.jpg"))

    rows = read_tsv(os.path.join(SI, "figS8_reference_validation.dat"))
    energy_rows = [r for r in rows if r["goal"] == "energy"]
    force_rows = [r for r in rows if r["goal"] == "force"]
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.75))
    x = np.arange(len(energy_rows))
    axes[0].bar(
        x,
        [fval(r["reference_vs_truth_error"]) for r in energy_rows],
        color=COLORS["blue"],
    )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([syslabel(r["system"]) for r in energy_rows])
    axes[0].set_ylabel(r"Reference to validation difference (meV $\cdot$ atom$^{-1}$)", fontsize=12)
    axes[0].set_title("Energy difference", fontsize=PANEL_TITLE_SIZE)
    style_axis(axes[0], grid="y")
    panel_label(axes[0], "a")

    x = np.arange(len(force_rows))
    axes[1].bar(
        x,
        [fval(r["reference_vs_truth_error"]) for r in force_rows],
        color=COLORS["orange"],
    )
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([syslabel(r["system"]) for r in force_rows])
    axes[1].set_ylabel(r"Reference to validation difference (eV $\cdot$ Å$^{-1}$)", fontsize=12)
    axes[1].set_title("Force difference", fontsize=PANEL_TITLE_SIZE)
    style_axis(axes[1], grid="y")
    panel_label(axes[1], "b")

    all_rows = energy_rows + force_rows
    x = np.arange(len(all_rows))
    ratios = [
        fval(r["truth_wall_time_s"]) / fval(r["reference_wall_time_s"])
        for r in all_rows
    ]
    colors = [COLORS["blue"]] * len(energy_rows) + [COLORS["orange"]] * len(force_rows)
    axes[2].bar(x, ratios, color=colors)
    axes[2].set_ylim(0, 4)
    axes[2].legend(
        handles=[
            Patch(facecolor=COLORS["blue"], label="Energy"),
            Patch(facecolor=COLORS["orange"], label="Force"),
        ],
        frameon=False,
        loc="upper left",
    )
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(
        [syslabel(r["system"]) for r in all_rows], rotation=25, ha="right"
    )
    axes[2].set_ylabel(
        r"Validation/reference wall-time ratio, $t_{\mathrm{val}}/t_{\mathrm{ref}}$", fontsize=12
    )
    axes[2].set_title("Validation cost", fontsize=PANEL_TITLE_SIZE)
    style_axis(axes[2], grid="y")
    panel_label(axes[2], "c")
    fig.subplots_adjust(wspace=0.46, bottom=0.20, top=0.88)
    savefig(fig, os.path.join(FIG_SI, "FigS2_reference_validation.jpg"))

    e = read_tsv(os.path.join(MAIN, "fig2_cost_to_target.dat"))
    f = read_tsv(os.path.join(SI, "figS3_force_cost.dat"))
    fig, ax = plt.subplots(figsize=(7.1, 4.8))
    for data, target_key, error_key, color, marker, label in [
        (
            e,
            "target_meV_per_atom",
            "actual_error_meV_per_atom",
            COLORS["blue"],
            "o",
            "Energy",
        ),
        (f, "target_eV_per_A", "actual_error_eV_per_A", COLORS["orange"], "s", "Force"),
    ]:
        x = [fval(r[error_key]) / fval(r[target_key]) for r in data]
        y = [fval(r["uniform_over_adaptive_speedup"]) for r in data]
        ax.scatter(x, y, s=46, color=color, marker=marker, label=label, zorder=3)
    ax.axvline(1.0, linestyle="--", color=COLORS["black"], linewidth=1.0)
    ax.axhline(1.0, linestyle="--", color=COLORS["black"], linewidth=1.0)
    ax.set_xlabel(r"Normalized realized error, $e_Q/\tau_Q$")
    ax.set_ylabel(r"Wall-time ratio, $t_{\mathrm{uniform}}/t_{\mathrm{GoalDFT}}$")
    ax.text(
        0.04,
        0.95,
        "target met",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=ANNOTATION_SIZE,
        color=COLORS["gray"],
    )
    ax.text(
        0.94,
        0.95,
        "GoalDFT faster",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=ANNOTATION_SIZE,
        color=COLORS["gray"],
    )
    ax.legend(frameon=False)
    style_axis(ax)
    savefig(fig, os.path.join(FIG_SI, "FigS3_accuracy_cost_map.jpg"))

    rows_out = []
    for r in e:
        rows_out.append(
            [
                syslabel(r["system"]),
                "Energy",
                r"meV $\cdot$ atom$^{-1}$",
                r["target_meV_per_atom"],
                fixed(fval(r["actual_error_meV_per_atom"]), 3),
                time_text(fval(r["adaptive_cost_s"])),
                time_text(fval(r["uniform_cost_s"])),
                fixed(fval(r["uniform_over_adaptive_speedup"]), 2),
                fixed(fval(r["final_Ecut_Ha"]), 1),
                str(ival(r["final_kgrid_n"])),
                sci(fval(r["final_scf_tol"]), 1),
            ]
        )
    for r in f:
        rows_out.append(
            [
                syslabel(r["system"]),
                "Force",
                r"eV $\cdot$ Å$^{-1}$",
                r["target_eV_per_A"],
                fixed(fval(r["actual_error_eV_per_A"]), 4),
                time_text(fval(r["adaptive_cost_s"])),
                time_text(fval(r["uniform_cost_s"])),
                fixed(fval(r["uniform_over_adaptive_speedup"]), 2),
                fixed(fval(r["final_Ecut_Ha"]), 1),
                str(ival(r["final_kgrid_n"])),
                sci(fval(r["final_scf_tol"]), 1),
            ]
        )
    tsv_cols = [
        "System",
        "Goal",
        "Unit",
        "Target",
        "Realized error",
        "GoalDFT wall time (s)",
        "Uniform wall time (s)",
        "Uniform/GoalDFT wall-time ratio",
        "Final Ecut (Ha)",
        "Final k-grid n",
        "Final SCF tolerance",
    ]
    display_cols = [
        "System",
        "Goal",
        "Unit",
        "Target",
        "Realized\nerror",
        "GoalDFT\nwall time (s)",
        "Uniform\nwall time (s)",
        r"$t_{\mathrm{uniform}}/t_{\mathrm{GoalDFT}}$",
        r"Final $E_{\mathrm{cut}}$" + "\n(Ha)",
        r"Final $n_k$",
        r"Final $\varepsilon_{\mathrm{SCF}}$",
    ]
    tsv_rows = []
    for row in rows_out:
        rr = list(row)
        rr[0] = rr[0].replace(r"TiO$_2$", "TiO2")
        rr[2] = (
            rr[2]
            .replace(r"meV $\cdot$ atom$^{-1}$", r"meV $\cdot$ atom^-1")
            .replace(r"eV $\cdot$ Å$^{-1}$", r"eV $\cdot$ Å^-1")
        )
        tsv_rows.append(rr)
    with open(
        os.path.join(TAB_SI, "TableS1_full_targets.tsv"),
        "w",
        encoding="utf-8",
        newline="",
    ) as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(tsv_cols)
        w.writerows(tsv_rows)
    write_table(
        tsv_cols,
        display_cols,
        rows_out,
        os.path.join(TAB_SI, "TableS1_full_targets_display.tsv"),
        os.path.join(TAB_SI, "TableS1_full_targets.jpg"),
        "Table S1. Complete GoalDFT target-by-target results.",
    )
    os.remove(os.path.join(TAB_SI, "TableS1_full_targets_display.tsv"))


if __name__ == "__main__":
    fig2_reliability()
    fig3_cost()
    fig4_calibration()
    fig5_budget_paths()
    table_main()
    figs_si()
