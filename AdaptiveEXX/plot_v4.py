import os, csv, math

os.chdir(os.path.split(os.path.realpath(__file__))[0])
import numpy as np
import matplotlib as mpl

mpl.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon

COLORS = {
    "blue": (31 / 255, 120 / 255, 180 / 255),
    "green": (51 / 255, 160 / 255, 44 / 255),
    "red": (227 / 255, 26 / 255, 28 / 255),
    "orange": (255 / 255, 127 / 255, 0 / 255),
    "purple": (106 / 255, 61 / 255, 154 / 255),
    "brown": (177 / 255, 89 / 255, 40 / 255),
    "gray": (0.38, 0.38, 0.38),
    "lightgray": (0.88, 0.88, 0.88),
    "black": (0.05, 0.05, 0.05),
}
COLOR_ORDER = [
    COLORS["blue"],
    COLORS["orange"],
    COLORS["green"],
    COLORS["purple"],
    COLORS["red"],
    COLORS["brown"],
    COLORS["gray"],
]
FONT_FAMILY = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
AXIS_LABEL_SIZE = 13.0
TICK_LABEL_SIZE = 11.0
LEGEND_SIZE = 10.2
PANEL_TITLE_SIZE = 12.0
PANEL_LABEL_SIZE = 12.5
ANNOTATION_SIZE = 10.0
FLOW_TEXT_SIZE = 10.6
FLOW_MATH_SIZE = 11.2
TABLE_TITLE_SIZE = 11.2
TABLE_FONT_SIZE = 9.0
LINE_WIDTH = 2.0
MARKER_SIZE = 6.0
SPINE_WIDTH = 0.9
DPI = 600
DEFECT_TOL = 3e-3

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

MAIN = os.path.join("data", "main")
SI = os.path.join("data", "si")
FIG_MAIN = os.path.join("figures", "main")
FIG_SI = os.path.join("figures", "si")
TAB_MAIN = os.path.join("tables", "main")
TAB_SI = os.path.join("tables", "si")
for d in [FIG_MAIN, FIG_SI, TAB_MAIN, TAB_SI]:
    os.makedirs(d, exist_ok=True)
    for name in os.listdir(d):
        if name.lower().endswith((".jpg", ".jpeg", ".png", ".tsv")):
            os.remove(os.path.join(d, name))

METHOD_ORDER = ["Vanilla_EXX", "ACE_occupied", "ACE_full", "AdaptiveEXX"]
METHOD_LABELS = {
    "Vanilla_EXX": "Direct EXX",
    "ACE_occupied": "ACE (occupied subspace)",
    "ACE_full": "ACE (extended subspace)",
    "AdaptiveEXX": "Adaptive EXX",
    "PeriodicACE_2": r"Periodic ACE ($p=2$)",
    "PeriodicACE_3": r"Periodic ACE ($p=3$)",
    "PeriodicACE_5": r"Periodic ACE ($p=5$)",
}

SYSTEM_LABELS = {
    "Si2": "Si (2 atoms)",
    "C2": "C (2 atoms)",
    "MgO2": "MgO (2 atoms)",
    "LiF2": "LiF (2 atoms)",
    "Si8": "Si (8 atoms)",
    "Si16": "Si (16 atoms)",
    "Si2_HSE06": "Si (2 atoms)",
    "MgO2_HSE06": "MgO (2 atoms)",
}
SYSTEM_TICKS = {
    "Si2": "Si (2)",
    "C2": "C (2)",
    "MgO2": "MgO (2)",
    "LiF2": "LiF (2)",
    "Si8": "Si (8)",
    "Si16": "Si (16)",
}


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


def energy_dev(r):
    return fval(
        r.get(
            "energy_deviation_from_vanilla_meV_per_atom",
            r.get("energy_error_meV_per_atom", "nan"),
        )
    )


def syslabel(name, short=False):
    return (SYSTEM_TICKS if short else SYSTEM_LABELS).get(
        name, name.replace("_HSE06", "")
    )


def savefig(fig, path):
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


def sci(x, digits=2):
    if not math.isfinite(x):
        return "NA"
    if x == 0:
        return "0"
    return f"{x:.{digits}e}"


def fixed(x, digits=2):
    if not math.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def time_text(x):
    if not math.isfinite(x):
        return "NA"
    return f"{x:.2f}" if x < 100 else f"{x:.1f}"


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
        -0.105,
        1.015,
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


def add_round_box(
    ax, xy, w, h, text, fc="white", ec=None, fontsize=FLOW_TEXT_SIZE, lw=1.0
):
    ec = COLORS["black"] if ec is None else ec
    p = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.018",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(p)
    ax.text(
        xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=fontsize
    )


def add_diamond(ax, center, w, h, text, fc="white", fontsize=FLOW_MATH_SIZE):
    x, y = center
    pts = np.array([[x, y + h / 2], [x + w / 2, y], [x, y - h / 2], [x - w / 2, y]])
    p = Polygon(
        pts, closed=True, facecolor=fc, edgecolor=COLORS["black"], linewidth=1.0
    )
    ax.add_patch(p)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize)


def add_arrow(ax, p1, p2, rad=0.0, text=None, text_xy=None):
    ax.add_patch(
        FancyArrowPatch(
            p1,
            p2,
            arrowstyle="-|>",
            mutation_scale=11,
            linewidth=1.0,
            color=COLORS["black"],
            connectionstyle=f"arc3,rad={rad}",
        )
    )
    if text is not None and text_xy is not None:
        ax.text(
            text_xy[0],
            text_xy[1],
            text,
            fontsize=ANNOTATION_SIZE,
            ha="center",
            va="center",
        )


def fig2_calibration():
    rows = [
        r
        for r in read_tsv(os.path.join(MAIN, "fig1_defect_calibration.dat"))
        if fval(r["full_occupied_defect"]) > 0
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.05), sharey=True)
    system_order = ["Si2", "MgO2", "Si16"]
    systems = [s for s in system_order if any(r["system"] == s for r in rows)]
    panels = [
        (
            "probe1_defect",
            None,
            r"Single-probe indicator, $\eta_1$",
            "Production indicator",
        ),
        (
            "probe3_defect",
            "probe3_stderr",
            r"Three-probe indicator, $\eta_3$",
            "Three-probe reference",
        ),
    ]
    for j, (ax, (field, errfield, xlabel, title)) in enumerate(zip(axes, panels)):
        vals = []
        for i, system in enumerate(systems):
            rr = [r for r in rows if r["system"] == system and fval(r[field]) > 0]
            x = np.asarray([fval(r[field]) for r in rr])
            y = np.asarray([fval(r["full_occupied_defect"]) for r in rr])
            vals.extend(x.tolist() + y.tolist())
            c = COLOR_ORDER[i]
            if errfield is None:
                ax.plot(
                    x,
                    y,
                    linestyle="none",
                    marker="o",
                    markersize=MARKER_SIZE,
                    color=c,
                    label=syslabel(system),
                )
            else:
                xe = np.asarray([2.0 * max(fval(r[errfield]), 0.0) for r in rr])
                ax.errorbar(
                    x,
                    y,
                    xerr=xe,
                    fmt="o",
                    markersize=MARKER_SIZE,
                    color=c,
                    ecolor=c,
                    elinewidth=0.9,
                    capsize=2,
                    label=syslabel(system),
                )
        if vals:
            lo, hi = min(vals), max(vals)
            ax.plot(
                [lo, hi],
                [lo, hi],
                "--",
                color=COLORS["black"],
                linewidth=1.0,
                label=r"$y=x$",
            )
        ax.axvline(DEFECT_TOL, linestyle=":", color=COLORS["gray"], linewidth=1.0)
        ax.axhline(DEFECT_TOL, linestyle=":", color=COLORS["gray"], linewidth=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=PANEL_TITLE_SIZE)
        style_axis(ax)
        panel_label(ax, chr(ord("a") + j))
    axes[0].set_ylabel(r"Full occupied-space defect, $\eta_{\mathrm{full}}$")
    axes[0].legend(frameon=False, handletextpad=0.45, labelspacing=0.35)
    axes[0].text(
        0.97,
        0.06,
        r"$\tau_{\mathrm{EXX}}=3\times10^{-3}$",
        transform=axes[0].transAxes,
        ha="right",
        va="bottom",
        fontsize=ANNOTATION_SIZE,
        color=COLORS["gray"],
    )
    axes[1].text(
        0.97,
        0.06,
        r"error bars: $\pm2\,\mathrm{SE}$",
        transform=axes[1].transAxes,
        ha="right",
        va="bottom",
        fontsize=ANNOTATION_SIZE,
        color=COLORS["gray"],
    )
    fig.subplots_adjust(wspace=0.18)
    savefig(fig, os.path.join(FIG_MAIN, "Fig2_defect_calibration.jpg"))


def fig3_scaling():
    rows = [
        r
        for r in read_tsv(os.path.join(MAIN, "fig4_size_scaling.dat"))
        if bval(r["converged"])
    ]
    atoms = sorted(set(ival(r["atoms"]) for r in rows))
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.25))
    colors = {"ACE_occupied": COLORS["blue"], "AdaptiveEXX": COLORS["orange"]}
    for method in ["ACE_occupied", "AdaptiveEXX"]:
        rr = sorted(
            [r for r in rows if r["method"] == method], key=lambda r: ival(r["atoms"])
        )
        axes[0].plot(
            [ival(r["atoms"]) for r in rr],
            [fval(r["total_wall_time_s"]) for r in rr],
            marker="o",
            markersize=MARKER_SIZE,
            linewidth=LINE_WIDTH,
            color=colors[method],
            label=METHOD_LABELS[method],
        )
    speed = []
    for n in atoms:
        ta = next(
            fval(r["total_wall_time_s"])
            for r in rows
            if ival(r["atoms"]) == n and r["method"] == "ACE_occupied"
        )
        td = next(
            fval(r["total_wall_time_s"])
            for r in rows
            if ival(r["atoms"]) == n and r["method"] == "AdaptiveEXX"
        )
        speed.append(ta / td)
    axes[1].plot(
        atoms,
        speed,
        marker="o",
        markersize=MARKER_SIZE,
        linewidth=LINE_WIDTH,
        color=COLORS["blue"],
    )
    axes[1].axhline(1.0, linestyle="--", color=COLORS["black"], linewidth=1.1)
    ymin = min(0.70, min(speed) - 0.08)
    ymax = max(speed) + 0.22
    axes[1].set_ylim(ymin, ymax)
    for n, s in zip(atoms, speed):
        dy = 10 if s >= 1 else -16
        axes[1].annotate(
            f"{s:.2f}×",
            (n, s),
            xytext=(0, dy),
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=ANNOTATION_SIZE,
            clip_on=False,
        )
    axes[0].set_xlabel(r"Number of atoms, $N$")
    axes[0].set_ylabel(r"Hybrid-stage wall time, $t_{\mathrm{hyb}}$ (s)")
    axes[0].set_yscale("log")
    axes[0].set_xticks(atoms)
    axes[0].legend(frameon=False, loc="upper left")
    axes[0].set_title("Wall-time crossover", fontsize=PANEL_TITLE_SIZE)
    axes[1].set_xlabel(r"Number of atoms, $N$")
    axes[1].set_ylabel(r"Speedup, $t_{\mathrm{ACE,occ}}/t_{\mathrm{AdaptiveEXX}}$")
    axes[1].set_xticks(atoms)
    axes[1].set_title("AdaptiveEXX speedup", fontsize=PANEL_TITLE_SIZE)
    for i, ax in enumerate(axes):
        style_axis(ax)
        panel_label(ax, chr(ord("a") + i))
    fig.subplots_adjust(wspace=0.30, left=0.10, right=0.98, bottom=0.16, top=0.90)
    savefig(fig, os.path.join(FIG_MAIN, "Fig3_size_crossover.jpg"))


def fig4_probe_efficiency():
    rows = read_tsv(os.path.join(MAIN, "fig4_probe_efficiency.dat"))
    order = ["Fixed 1", "Fixed 3", "Fixed 5"]
    systems = [s for s in ["Si8", "MgO2"] if any(r.get("system") == s for r in rows)]
    fig, axes = plt.subplots(1, 3, figsize=(11.8, 4.15))
    x = np.arange(len(order))
    width = 0.34 if len(systems) > 1 else 0.58
    system_colors = [COLORS["blue"], COLORS["orange"]]
    for i, system in enumerate(systems):
        rr0 = [r for r in rows if r.get("system") == system]
        rr = [next((r for r in rr0 if r["strategy"] == name), None) for name in order]
        fixed1 = next(
            (fval(r["total_wall_time_s"]) for r in rr0 if r["strategy"] == "Fixed 1"),
            np.nan,
        )
        rel_time = [
            (
                fval(r["total_wall_time_s"]) / fixed1
                if r is not None and fixed1 > 0
                else np.nan
            )
            for r in rr
        ]
        probes = [
            fval(r["total_monitoring_probes"]) if r is not None else np.nan for r in rr
        ]
        evals = [fval(r["scf_evaluations"]) if r is not None else np.nan for r in rr]
        offset = (i - (len(systems) - 1) / 2) * width
        axes[0].bar(
            x + offset,
            rel_time,
            width=width,
            color=system_colors[i],
            label=syslabel(system),
        )
        axes[1].bar(x + offset, probes, width=width, color=system_colors[i])
        axes[2].bar(x + offset, evals, width=width, color=system_colors[i])
    axes[0].axhline(1.0, linestyle="--", color=COLORS["black"], linewidth=1.0)
    axes[0].set_ylabel(r"Relative wall time, $t_{n_p}/t_1$")
    axes[1].set_ylabel("Total stochastic probes")
    axes[2].set_ylabel("SCF evaluations")
    titles = ["Monitoring overhead", "Stochastic work", "SCF effort"]
    for j, ax in enumerate(axes):
        ax.set_xticks(x)
        ax.set_xticklabels(["1", "3", "5"])
        ax.set_xlabel(r"Probes per monitored call, $n_p$")
        ax.set_title(titles[j], fontsize=PANEL_TITLE_SIZE)
        style_axis(ax, grid="y")
        panel_label(ax, chr(ord("a") + j))
    if systems:
        axes[0].legend(frameon=False)
    fig.subplots_adjust(wspace=0.34, left=0.075, right=0.985, bottom=0.18, top=0.89)
    savefig(fig, os.path.join(FIG_MAIN, "Fig4_probe_count_sensitivity.jpg"))


def fig5_controller():
    summary = [
        r
        for r in read_tsv(os.path.join(MAIN, "fig5_controller_efficiency.dat"))
        if r["functional"] == "PBE0"
    ]
    order = ["Si2", "C2", "MgO2", "LiF2", "Si8"]
    summary = [
        next(r for r in summary if r["system"] == s)
        for s in order
        if any(x["system"] == s for x in summary)
    ]
    systems = [r["system"] for r in summary]
    x = np.arange(len(summary))

    trace = [
        r
        for r in read_tsv(os.path.join(SI, "figS10_adaptive_trace.dat"))
        if r.get("functional") == "PBE0" and r.get("system") in systems
    ]

    initial, cap, defect = [], [], []
    for system in systems:
        rr = sorted(
            [r for r in trace if r["system"] == system],
            key=lambda r: ival(r["call"]),
        )
        n_initial = sum(
            1 for r in rr if r["action"] == "rebuild" and ival(r["call"]) == 1
        )
        n_defect = sum(
            1
            for r in rr
            if r["action"] == "rebuild"
            and ival(r["call"]) != 1
            and fval(r.get("probes_used", "0")) > 0
        )
        n_cap = sum(
            1
            for r in rr
            if r["action"] == "rebuild"
            and ival(r["call"]) != 1
            and fval(r.get("probes_used", "0")) <= 0
        )
        initial.append(n_initial)
        cap.append(n_cap)
        defect.append(n_defect)

    reuses = np.asarray([fval(r["reuses"]) for r in summary])
    initial = np.asarray(initial, dtype=float)
    cap = np.asarray(cap, dtype=float)
    defect = np.asarray(defect, dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.25))

    axes[0].bar(x, initial, color=COLORS["gray"], label="Initial build")
    axes[0].bar(
        x, cap, bottom=initial, color=COLORS["orange"], label="Safety-cap rebuild"
    )
    axes[0].bar(
        x,
        defect,
        bottom=initial + cap,
        color=COLORS["red"],
        label="Defect-triggered rebuild",
    )
    axes[0].bar(
        x,
        reuses,
        bottom=initial + cap + defect,
        color=COLORS["green"],
        label="Reuse",
    )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([syslabel(s, short=True) for s in systems])
    axes[0].set_ylabel("Controller actions")
    axes[0].set_title("Controller action statistics", fontsize=PANEL_TITLE_SIZE)

    totals = initial + cap + defect + reuses
    axes[0].set_ylim(0, max(36.0, float(np.max(totals)) * 1.18))

    axes[0].legend(
        frameon=False,
        ncol=2,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.00),
        columnspacing=1.0,
        handletextpad=0.45,
    )
    style_axis(axes[0], grid="y")
    panel_label(axes[0], "a")

    probes = [fval(r["total_monitoring_probes"]) for r in summary]
    exx_time = [fval(r["total_exx_call_time_s"]) for r in summary]
    bars = axes[1].bar(
        x, probes, width=0.58, color=COLORS["blue"], label="Stochastic probes"
    )
    ax2 = axes[1].twinx()
    (line,) = ax2.plot(
        x,
        exx_time,
        marker="o",
        markersize=MARKER_SIZE,
        linewidth=LINE_WIDTH,
        color=COLORS["orange"],
        label="EXX call time",
    )
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([syslabel(s, short=True) for s in systems])
    axes[1].set_ylabel("Total stochastic probes", color=COLORS["blue"])
    ax2.set_ylabel("Total EXX call time (s)", color=COLORS["orange"])
    ax2.set_yscale("log")
    style_twin_axis(ax2, COLORS["orange"])
    axes[1].set_title("Monitoring and EXX cost", fontsize=PANEL_TITLE_SIZE)
    style_axis(axes[1], grid="y")
    panel_label(axes[1], "b")
    axes[1].legend(
        [bars, line],
        ["Stochastic probes", "EXX call time"],
        frameon=False,
        loc="upper left",
    )

    fig.subplots_adjust(wspace=0.36, left=0.085, right=0.91, bottom=0.18, top=0.90)
    savefig(fig, os.path.join(FIG_MAIN, "Fig5_controller_statistics.jpg"))


def table_main():
    perf = read_tsv(os.path.join(MAIN, "table2_performance.dat"))
    ctrl = {
        r["system"]: r
        for r in read_tsv(os.path.join(MAIN, "table3_controller_statistics.dat"))
    }
    nat = {"Si2": 2, "C2": 2, "MgO2": 2, "LiF2": 2, "Si8": 8}
    rows = []
    for system in ["Si2", "C2", "MgO2", "LiF2", "Si8"]:
        direct = next(
            r for r in perf if r["system"] == system and r["method"] == "Vanilla_EXX"
        )
        ace = next(
            r for r in perf if r["system"] == system and r["method"] == "ACE_occupied"
        )
        adp = next(
            r for r in perf if r["system"] == system and r["method"] == "AdaptiveEXX"
        )
        c = ctrl[system]
        dev = (
            abs(fval(adp["final_energy_Ha"]) - fval(direct["final_energy_Ha"]))
            * 27211.386245988
            / nat[system]
        )
        rows.append(
            [
                syslabel(system),
                time_text(fval(adp["total_wall_time_s"])),
                time_text(fval(ace["total_wall_time_s"])),
                fixed(
                    fval(ace["total_wall_time_s"]) / fval(adp["total_wall_time_s"]), 2
                ),
                sci(dev, 2),
                str(ival(c["rebuilds"])),
                str(ival(c["reuses"])),
                str(ival(c["total_monitoring_probes"])),
            ]
        )
    tsv_cols = [
        "System",
        "AdaptiveEXX wall time (s)",
        "Occupied ACE wall time (s)",
        "ACE/AdaptiveEXX wall-time ratio",
        "Energy deviation from direct EXX (meV atom^-1)",
        "Rebuilds",
        "Reuses",
        "Monitoring probes",
    ]
    display_cols = [
        "System",
        "AdaptiveEXX\nwall time (s)",
        "ACE (occupied)\nwall time (s)",
        "Wall-time ratio\n" + r"$t_{\mathrm{ACE,occ}}/t_{\mathrm{AdaptiveEXX}}$",
        r"$|\Delta E|$ vs direct EXX" + "\n" + r"(meV atom$^{-1}$)",
        "Rebuilds",
        "Reuses",
        "Probes",
    ]
    write_table(
        tsv_cols,
        display_cols,
        rows,
        os.path.join(TAB_MAIN, "Table1_summary.tsv"),
        os.path.join(TAB_MAIN, "Table1_summary.jpg"),
        "Table 1. PBE0 performance and controller statistics.",
    )


def figs_si():
    rows = read_tsv(os.path.join(SI, "figS4_all_pbe_convergence.dat"))
    order = ["Si2", "C2", "MgO2", "LiF2", "Si8"]
    systems = [s for s in order if any(r["system"] == s for r in rows)]
    ncols = 3
    nrows = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.6, 6.9), sharey=True)
    axes = np.array(axes).reshape(-1)
    for idx, (ax, system) in enumerate(zip(axes, systems)):
        for i, method in enumerate(METHOD_ORDER):
            rr = [r for r in rows if r["system"] == system and r["method"] == method]
            rr.sort(key=lambda r: ival(r["evaluation"]))
            if rr:
                ax.plot(
                    [ival(r["evaluation"]) for r in rr],
                    [fval(r["density_residual"]) for r in rr],
                    linewidth=1.9,
                    color=COLOR_ORDER[i],
                    label=METHOD_LABELS[method],
                    marker="o",
                    markersize=3.2,
                    markevery=max(1, len(rr) // 8),
                )
        ax.set_yscale("log")
        ax.set_title(syslabel(system), fontsize=PANEL_TITLE_SIZE)
        if idx // ncols == nrows - 1:
            ax.set_xlabel("SCF evaluation")
        if idx % ncols == 0:
            ax.set_ylabel(r"Density residual, $\|\Delta\rho\|$")
        style_axis(ax)

    leg_ax = axes[-1]
    leg_ax.axis("off")
    h, l = axes[0].get_legend_handles_labels()
    leg_ax.legend(
        h,
        l,
        ncol=1,
        frameon=False,
        loc="center",
        fontsize=LEGEND_SIZE + 0.5,
        handlelength=2.6,
        labelspacing=0.75,
    )
    # leg_ax.text(
    # 0.5,
    # 0.77,
    # "Methods",
    # transform=leg_ax.transAxes,
    # ha="center",
    # va="center",
    # fontsize=PANEL_TITLE_SIZE,
    # fontweight="normal",
    # )
    fig.subplots_adjust(
        hspace=0.34, wspace=0.25, left=0.075, right=0.985, bottom=0.10, top=0.94
    )
    savefig(fig, os.path.join(FIG_SI, "FigS1_all_pbe_convergence.jpg"))

    rows = [
        r
        for r in read_tsv(os.path.join(SI, "figS2_sensitivity.dat"))
        if r.get("system") == "Si8"
    ]

    def hit(tau, npb):
        return next(
            r
            for r in rows
            if abs(fval(r["defect_tolerance"]) - tau) < 1e-12
            and ival(r["probe_count"]) == npb
        )

    display_rows = [
        hit(1e-2, 3),
        hit(3e-3, 3),
        hit(1e-3, 3),
        hit(3e-3, 1),
        hit(3e-3, 3),
        hit(3e-3, 5),
    ]
    x = np.asarray([0, 1, 2, 4, 5, 6], dtype=float)
    labels = [r"$10^{-2}$", r"$3\times10^{-3}$", r"$10^{-3}$", "1", "3", "5"]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.35))
    ysets = [
        (
            [fval(r["total_wall_time_s"]) for r in display_rows],
            r"Hybrid-stage wall time, $t_{\mathrm{hyb}}$ (s)",
            COLORS["blue"],
            "Wall-time sensitivity",
        ),
        (
            [fval(r["rebuilds"]) for r in display_rows],
            "ACE rebuilds",
            COLORS["orange"],
            "Rebuild sensitivity",
        ),
    ]
    for j, (ax, (vals, ylabel, color, title)) in enumerate(zip(axes, ysets)):
        ax.bar(x, vals, width=0.72, color=color)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=0)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=PANEL_TITLE_SIZE)
        style_axis(ax, grid="y")
        panel_label(ax, chr(ord("a") + j))
        ax.axvline(3.0, color=COLORS["lightgray"], linewidth=1.0)
        group_label_y = -0.235
        ax.text(
            1.0,
            group_label_y,
            r"Threshold ($n_p=3$)",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="baseline",
            fontsize=TICK_LABEL_SIZE,
            clip_on=False,
        )
        ax.text(
            5.0,
            group_label_y,
            r"Probe count ($\tau_{\mathrm{EXX}}=3\times10^{-3}$)",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="baseline",
            fontsize=TICK_LABEL_SIZE,
            clip_on=False,
        )
    fig.subplots_adjust(bottom=0.30, wspace=0.30, left=0.09, right=0.985, top=0.90)
    savefig(fig, os.path.join(FIG_SI, "FigS2_sensitivity.jpg"))

    rows = read_tsv(os.path.join(SI, "figS7_calibration_relative_error.dat"))
    order = ["Si2", "MgO2", "Si16"]
    grouped = [s for s in order if any(r["system"] == s for r in rows)]
    seed_data = [
        [
            fval(r["independent_probe1_relative_error"])
            for r in rows
            if r["system"] == s and fval(r["independent_probe1_relative_error"]) > 0
        ]
        for s in grouped
    ]
    med1, med3 = [], []
    for system in grouped:
        rr = [r for r in rows if r["system"] == system]
        e1 = [
            fval(r["probe1_relative_error"])
            for r in rr
            if fval(r["probe1_relative_error"]) > 0
        ]
        e3 = [
            fval(r["probe3_relative_error"])
            for r in rr
            if fval(r["probe3_relative_error"]) > 0
        ]
        med1.append(float(np.median(e1)) if e1 else np.nan)
        med3.append(float(np.median(e3)) if e3 else np.nan)
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.25))
    bp = axes[0].boxplot(
        seed_data,
        patch_artist=True,
        tick_labels=[syslabel(s, short=True) for s in grouped],
        showfliers=False,
    )
    for i, p in enumerate(bp["boxes"]):
        p.set_facecolor(COLOR_ORDER[i])
        p.set_edgecolor(COLORS["black"])
    for key in ["whiskers", "caps", "medians"]:
        for a in bp[key]:
            a.set_color(COLORS["black"])
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Independent-seed relative error")
    axes[0].set_title("Single-probe seed robustness", fontsize=PANEL_TITLE_SIZE)
    style_axis(axes[0], grid="y")
    panel_label(axes[0], "a")

    x = np.arange(len(grouped))
    width = 0.34
    axes[1].bar(x - width / 2, med1, width=width, color=COLORS["blue"], label="1 probe")
    axes[1].bar(
        x + width / 2, med3, width=width, color=COLORS["orange"], label="3 probes"
    )
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([syslabel(s, short=True) for s in grouped])
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Median relative error")
    axes[1].set_title("Probe-count accuracy", fontsize=PANEL_TITLE_SIZE)
    axes[1].legend(frameon=False, loc="upper right")
    style_axis(axes[1], grid="y")
    panel_label(axes[1], "b")
    fig.subplots_adjust(wspace=0.30, left=0.095, right=0.985, bottom=0.16, top=0.90)
    savefig(fig, os.path.join(FIG_SI, "FigS3_single_probe_robustness.jpg"))

    rows = read_tsv(os.path.join(MAIN, "fig6_hse_transfer.dat"))
    systems = [
        s for s in ["Si2_HSE06", "MgO2_HSE06"] if any(r["system"] == s for r in rows)
    ]
    fig, axes = plt.subplots(1, len(systems), figsize=(9.6, 4.15), sharey=True)
    axes = [axes] if len(systems) == 1 else axes
    for j, (ax, system) in enumerate(zip(axes, systems)):
        rr = [r for r in rows if r["system"] == system]
        for i, method in enumerate(METHOD_ORDER):
            hit0 = [r for r in rr if r["method"] == method]
            if hit0 and energy_dev(hit0[0]) > 0:
                ax.scatter(
                    fval(hit0[0]["inner_scf_time_s"]),
                    energy_dev(hit0[0]),
                    s=62,
                    color=COLOR_ORDER[i],
                    label=METHOD_LABELS[method],
                    zorder=3,
                )
        ax.set_yscale("log")
        ax.set_xlabel("Hybrid-SCF time (s)")
        ax.set_title(syslabel(system), fontsize=PANEL_TITLE_SIZE)
        style_axis(ax)
        panel_label(ax, chr(ord("a") + j))
        ax.margins(x=0.18)
    axes[0].set_ylabel(r"Energy deviation from direct EXX (meV atom$^{-1}$)")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(
        h, l, ncol=3, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.995)
    )
    fig.subplots_adjust(wspace=0.24, left=0.10, right=0.985, bottom=0.16, top=0.83)
    savefig(fig, os.path.join(FIG_SI, "FigS4_hse_transferability.jpg"))

    rows = read_tsv(os.path.join(SI, "figS9_periodic_rebuild.dat"))
    systems = [s for s in ["MgO2", "Si8"] if any(r["system"] == s for r in rows)]
    fig, axes = plt.subplots(1, len(systems), figsize=(10.8, 4.45))
    axes = [axes] if len(systems) == 1 else axes
    order = [
        "ACE_occupied",
        "AdaptiveEXX",
        "PeriodicACE_2",
        "PeriodicACE_3",
        "PeriodicACE_5",
    ]
    short_labels = {
        "ACE_occupied": "Occ.\nACE",
        "AdaptiveEXX": "Adaptive\nEXX",
        "PeriodicACE_2": "Periodic\n$p=2$",
        "PeriodicACE_3": "Periodic\n$p=3$",
        "PeriodicACE_5": "Periodic\n$p=5$",
    }
    for j, (ax, system) in enumerate(zip(axes, systems)):
        rr0 = [r for r in rows if r["system"] == system]
        rr = [
            next(r for r in rr0 if r["method"] == m)
            for m in order
            if any(x["method"] == m for x in rr0)
        ]
        x = np.arange(len(rr))
        vals = [fval(r["total_wall_time_s"]) for r in rr]
        ax.bar(x, vals, color=COLOR_ORDER[: len(rr)], width=0.68)
        ax.set_xticks(x)
        ax.set_xticklabels(
            [short_labels[r["method"]] for r in rr], rotation=0, fontsize=10.0
        )
        ax.set_ylabel(r"Hybrid-stage wall time, $t_{\mathrm{hyb}}$ (s)")
        ax.set_title(syslabel(system), fontsize=PANEL_TITLE_SIZE)
        style_axis(ax, grid="y")
        panel_label(ax, chr(ord("a") + j))
        ax.margins(x=0.06)
    fig.subplots_adjust(bottom=0.24, wspace=0.25, left=0.08, right=0.985, top=0.90)
    savefig(fig, os.path.join(FIG_SI, "FigS5_fixed_period_baselines.jpg"))

    perf = read_tsv(os.path.join(MAIN, "table2_performance.dat"))
    hse = read_tsv(os.path.join(SI, "tableS1_hse_performance.dat"))
    rows_out = []
    for r in perf + hse:
        rows_out.append(
            [
                syslabel(r["system"]),
                r["functional"],
                METHOD_LABELS.get(r["method"], r["method"]),
                "Yes" if bval(r["converged"]) else "No",
                str(ival(r["scf_evaluations"])),
                time_text(fval(r["total_wall_time_s"])),
                sci(fval(r["final_density_residual"]), 2),
            ]
        )
    tsv_cols = [
        "System",
        "Functional / case",
        "Method",
        "Converged",
        "SCF evaluations",
        "Wall time (s)",
        "Final density residual",
    ]
    display_cols = [
        "System",
        "Functional / case",
        "Method",
        "Converged",
        "SCF evals",
        "Wall time (s)",
        r"Final $\|\Delta\rho\|$",
    ]
    write_table(
        tsv_cols,
        display_cols,
        rows_out,
        os.path.join(TAB_SI, "TableS1_full_benchmark.tsv"),
        os.path.join(TAB_SI, "TableS1_full_benchmark.jpg"),
        "Table S1. Complete PBE0 and HSE06 benchmarks.",
    )


if __name__ == "__main__":
    fig2_calibration()
    fig3_scaling()
    fig4_probe_efficiency()
    fig5_controller()
    table_main()
    figs_si()
