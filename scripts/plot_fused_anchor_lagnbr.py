"""Draw the MERGED_OLS_LAGNBR figures and the report tables.

Reads the CSVs written by `scripts/run_fused_anchor_lagnbr_report.jl` and writes PNG (300 dpi) and PDF figures
to `output/fused_anchor_lagnbr/figures/`, then copies the figures and a set of rounded report tables into
`Interim_results/` so the Typst report compiles from that folder alone:

    py -3.13 scripts/plot_fused_anchor_lagnbr.py

Same visual system as `plot_gauge_satellite_fusion.py`. The two merged anchors are the only categorical pair:
MERGED_OLS in muted grey (the incumbent) and MERGED_OLS_LAGNBR in blue (the new one), so every figure reads
"blue against grey". Improvements use the blue (better) / red (worse) diverging ramp with a grey centre.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "output" / "fused_anchor_lagnbr"
FIG_DIR = RUN_DIR / "figures"
REPORT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "Interim_results"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
ORDINAL_BLUES = ["#86b6ef", "#5598e7", "#256abf", "#104281"]
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # validated categorical slots 1-4
NEW_COLOR = ORDINAL_BLUES[2]
OLD_COLOR = MUTED
WORSE_COLOR = "#e34948"
DIVERGING = LinearSegmentedColormap.from_list(
    "better_worse", ["#b52f2f", "#e34948", "#ef8d8c", "#f0efec", "#86b6ef", "#3987e5", "#1c5cab"])

NEW, OLD = "MERGED_OLS_LAGNBR", "MERGED_OLS"
ANCHORS = ["FY4B", "GPM", "GSMaP", "MERGED_MEAN", OLD, NEW]
ANCHOR_LABELS = {"FY4B": "FY4B", "GPM": "GPM", "GSMaP": "GSMaP", "MERGED_MEAN": "Merged\nequal",
                 OLD: "Merged\nOLS", NEW: "Merged\nLAGNBR"}
SCHEMES = [("balanced_spatial", "Balanced spatial folds"), ("random", "Random folds")]
METHOD_LABELS = {
    "raw": "Raw anchor", "adw": "ADW (gauges only)", "idw": "IDW", "tps": "TPS", "gwr": "GWR (gauges only)",
    "auto": "GWR family, chosen in-fold", "residual_gwr": "Residual GWR", "mixed_gwr": "Mixed GWR", "mgwr": "MGWR",
    "blend_residual_gwr": "Residual GWR + ADW blend", "blend_mixed_gwr": "Mixed GWR + ADW blend",
    "blend_mgwr": "MGWR + ADW blend", "blend_agrenv_residual_gwr": "Residual GWR + agreement blend",
    "blend_agrenv_mixed_gwr": "Mixed GWR + agreement blend", "blend_agrenv_mgwr": "MGWR + agreement blend",
}
# The methods that read the anchor. The gauge-only ones score identically on both anchors by construction.
ANCHORED = ["raw", "residual_gwr", "mixed_gwr", "mgwr", "auto", "blend_residual_gwr", "blend_mixed_gwr",
            "blend_mgwr", "blend_agrenv_residual_gwr", "blend_agrenv_mixed_gwr", "blend_agrenv_mgwr"]
ALL_METHODS = ["raw", "idw", "adw", "tps", "gwr"] + ANCHORED[1:]
VS_ADW_METHODS = ["residual_gwr", "mgwr", "auto", "blend_mgwr", "blend_agrenv_mgwr"]
STRATA = [("all", "All hours"), ("no_rain", "Dry\n< 0.1"), ("light", "Light\n0.1–2.5"),
          ("moderate", "Moderate\n2.5–8"), ("heavy", "Heavy\n≥ 8")]
STRATUM_NAMES = {"all": "all", "no_rain": "dry", "light": "light", "moderate": "moderate", "heavy": "heavy"}
TERMS = ["intercept", "fy4b", "gpm", "gsmap", "fy4b_lag1", "gpm_lag1", "gsmap_lag1", "fy4b_lead1", "gpm_lead1",
         "gsmap_lead1", "fy4b_nbr8", "gpm_nbr8", "gsmap_nbr8"]
ENVELOPE = ["0.1–2.5", "2.5–8", "≥ 8"]
TEXT_COLUMNS = ("scheme", "product", "method", "stratum", "baseline_product", "metric", "term", "file",
                "comparison", "status", "covariate_variables", "best_traditional")

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK_SECONDARY,
    "axes.titlesize": 7.5,
    "axes.titlecolor": INK,
    "text.color": INK,
    "xtick.color": AXIS,
    "ytick.color": AXIS,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "pdf.fonttype": 42,
})


def read_rows(path: Path) -> list[dict]:
    """CSV rows with numeric-looking fields converted to float and true/false to bool."""
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            parsed = {}
            for key, value in row.items():
                if value in ("true", "false"):
                    parsed[key] = value == "true"
                elif key in TEXT_COLUMNS:
                    parsed[key] = value
                else:
                    try:
                        parsed[key] = float(value)
                    except ValueError:
                        parsed[key] = value
            rows.append(parsed)
    return rows


def one(rows: list[dict], **conditions) -> dict | None:
    for row in rows:
        if all(row.get(k) == v for k, v in conditions.items()):
            return row
    return None


def save(fig: plt.Figure, name: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300)
    fig.savefig(FIG_DIR / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


def clean_axis(ax: plt.Axes, grid_axis: str = "y") -> None:
    ax.grid(axis=grid_axis, color=GRID, linewidth=0.5, zorder=0)
    ax.tick_params(length=2, pad=1.5, labelsize=6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def title(fig: plt.Figure, heading: str, subheading: str, top: float = 0.975) -> None:
    fig.suptitle(heading, x=0.02, ha="left", fontsize=9, y=top, va="top")
    fig.text(0.02, top - 0.24 / fig.get_figheight(), subheading, ha="left", va="top", fontsize=6.3,
             color=INK_SECONDARY, linespacing=1.3)


def pct(x: float) -> float:
    return 100.0 * x


def zero_line(ax: plt.Axes, vertical: bool = True) -> None:
    (ax.axvline if vertical else ax.axhline)(0, color=INK_SECONDARY, linewidth=0.6, zorder=1)


def ci_dot(ax, x, y, low, high, color, vertical=False, hollow=False, size=4.5):
    """A point with its 95% interval: a thin whisker and a marker ringed in the surface colour."""
    if vertical:
        ax.plot([x, x], [low, high], color=color, linewidth=1.1, zorder=2, solid_capstyle="round")
    else:
        ax.plot([low, high], [y, y], color=color, linewidth=1.1, zorder=2, solid_capstyle="round")
    ax.plot(x, y, marker="o", markersize=size, markerfacecolor=SURFACE if hollow else color,
            markeredgecolor=color if hollow else SURFACE, markeredgewidth=0.9 if hollow else 0.7, zorder=3)


# --------------------------------------------------------------------------------------------------
# Figure 1: the new anchor against MERGED_OLS, method by method and stratum by stratum


GAIN_LIMIT = 8.0   # % - the diverging ramp saturates here; the raw anchor's dry-hour loss runs past it


def gain_cell(ax, x, y, row, limit=GAIN_LIMIT):
    """One heat-strip cell: colour = gain, text = value, bold when the interval excludes zero."""
    gain = pct(row["relative_improvement"])
    significant = row["ci_low"] > 0 or row["ci_high"] < 0
    colour = DIVERGING(TwoSlopeNorm(0, -limit, limit)(np.clip(gain, -limit, limit)))
    ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor=colour, edgecolor=SURFACE, linewidth=1.2))
    shade = abs(gain) / limit > 0.55
    ax.text(x, y, f"{gain:+.1f}", ha="center", va="center", fontsize=4.9,
            fontweight="bold" if significant else "normal",
            color=SURFACE if shade else (INK if significant else MUTED))


def fig1_gain(comparison: list[dict]) -> None:
    fig = plt.figure(figsize=(7.0, 3.9))
    grid = fig.add_gridspec(1, 4, width_ratios=[1.2, 1.25, 1.2, 1.25], left=0.2, right=0.995, top=0.84,
                            bottom=0.13, wspace=0.08)
    title(fig, "The new anchor improves nearly every correction method, mostly in moderate rain",
          "RMSE improvement of MERGED_OLS_LAGNBR over MERGED_OLS (%), same method on both sides · whisker / bold = "
          "95% day-block interval excludes zero")
    positions = np.arange(len(ANCHORED))[::-1]
    for column, (scheme, heading) in enumerate(SCHEMES):
        ax = fig.add_subplot(grid[0, 2 * column])
        clean_axis(ax, "x")
        zero_line(ax)
        for y, method in zip(positions, ANCHORED):
            row = one(comparison, scheme=scheme, method=method, stratum="all")
            better, worse = row["ci_low"] > 0, row["ci_high"] < 0
            ci_dot(ax, pct(row["relative_improvement"]), y, pct(row["ci_low"]), pct(row["ci_high"]),
                   NEW_COLOR if better else WORSE_COLOR if worse else MUTED, hollow=not (better or worse), size=4.0)
        ax.set_xlim(-0.8, 2.3)
        ax.set_ylim(-0.6, len(ANCHORED) - 0.4)
        ax.set_title(heading.replace(" folds", ""), loc="left")
        ax.set_xlabel("All hours: RMSE improvement (%)")
        if column == 0:
            ax.set_yticks(positions, [METHOD_LABELS[m] for m in ANCHORED])
        else:
            ax.set_yticks(positions, [])

        strip = fig.add_subplot(grid[0, 2 * column + 1])
        for y, method in zip(positions, ANCHORED):
            row = one(comparison, scheme=scheme, method=method, stratum="all")
            strip.text(-0.95, y, f"{row['RMSE_baseline']:.3f}\n→{row['RMSE_treatment']:.3f}", ha="center",
                       va="center", fontsize=4.4, color=INK_SECONDARY, linespacing=1.1)
        for x, (stratum, _) in enumerate(STRATA[1:]):
            for y, method in zip(positions, ANCHORED):
                gain_cell(strip, x, y, one(comparison, scheme=scheme, method=method, stratum=stratum))
        strip.set_xlim(-1.6, len(STRATA) - 1.5)
        strip.set_ylim(-0.6, len(ANCHORED) - 0.4)
        strip.set_xticks([-0.95] + list(range(len(STRATA) - 1)),
                         ["RMSE"] + [label for _, label in STRATA[1:]], fontsize=4.7)
        strip.set_yticks([])
        strip.tick_params(length=0)
        for spine in strip.spines.values():
            spine.set_visible(False)
        strip.set_title("by gauge intensity, mm/h", loc="right", fontsize=6.0, color=INK_SECONDARY)
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=NEW_COLOR,
                      markeredgecolor=SURFACE, label="Better"),
               Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=WORSE_COLOR,
                      markeredgecolor=SURFACE, label="Worse"),
               Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=SURFACE,
                      markeredgecolor=MUTED, label="Interval includes zero"),
               Line2D([], [], linestyle="none", label="RMSE column: on OLS → on LAGNBR (mm/h)")]
    fig.legend(handles=handles, frameon=False, fontsize=6.2, ncol=4, loc="lower left", bbox_to_anchor=(0.2, 0.0))
    save(fig, "lagnbr_fig1_gain_by_method")


# --------------------------------------------------------------------------------------------------
# Figure 2: every anchor, on the cells all six share


SHORT_ANCHORS = {"FY4B": "FY4B", "GPM": "GPM", "GSMaP": "GSMaP", "MERGED_MEAN": "Equal", OLD: "OLS", NEW: "LAGNBR"}
PANEL = [("raw", "hollow", MUTED), ("residual_gwr", "filled", ORDINAL_BLUES[0]),
         ("mgwr", "filled", ORDINAL_BLUES[1]), ("blend_mgwr", "filled", ORDINAL_BLUES[2]),
         ("blend_agrenv_mgwr", "filled", ORDINAL_BLUES[3])]


def fig2_anchors(metrics: list[dict]) -> None:
    rows = [r for r in metrics if r["scheme"] == "balanced_spatial"]
    fig, axes = plt.subplots(2, len(STRATA), figsize=(7.0, 4.6))
    fig.subplots_adjust(left=0.07, right=0.995, top=0.85, bottom=0.15, wspace=0.34, hspace=0.12)
    title(fig, "The new merge is the best anchor overall and mostly in moderate rain; heavy rain favours the equal merge",
          "Six anchors on the 3.02 M shared cells, balanced spatial folds · ring = lowest RMSE for that method · "
          "shaded = new anchor · Equal, OLS = earlier merges")
    x = np.arange(len(ANCHORS))
    offsets = np.linspace(-0.3, 0.3, len(PANEL))
    for column, (stratum, heading) in enumerate(STRATA):
        for row_index, metric in enumerate(("RMSE", "Bias")):
            ax = axes[row_index, column]
            clean_axis(ax)
            ax.axvspan(len(ANCHORS) - 1.5, len(ANCHORS) - 0.5, color=GRID, alpha=0.55, zorder=0, linewidth=0)
            adw = one(rows, product="gauge_only", method="adw", stratum=stratum)[metric]
            ax.axhline(adw, color=INK, linestyle=(0, (3, 2)), linewidth=0.8, zorder=1)
            if metric == "Bias":
                zero_line(ax, vertical=False)
            for offset, (method, style, color) in zip(offsets, PANEL):
                values = [one(rows, product=a, method=method, stratum=stratum)[metric] for a in ANCHORS]
                ax.plot(x + offset, values, linestyle="none", marker="o", markersize=3.4,
                        markerfacecolor=SURFACE if style == "hollow" else color,
                        markeredgecolor=color if style == "hollow" else SURFACE,
                        markeredgewidth=0.8 if style == "hollow" else 0.5, zorder=3)
                if metric == "RMSE" and method != "raw":
                    best = int(np.argmin(values))
                    ax.plot(best + offset, values[best], marker="o", markersize=5.4, markerfacecolor="none",
                            markeredgecolor=INK, markeredgewidth=0.7, zorder=4)
            ax.set_xticks(x, [SHORT_ANCHORS[a] for a in ANCHORS] if row_index == 1 else [], fontsize=5.0,
                          rotation=60, ha="right", rotation_mode="anchor")
            ax.tick_params(labelsize=5.2)
            if row_index == 0:
                ax.set_title(heading.replace("\n", " "), loc="left", fontsize=6.8)
        if stratum == "all":
            axes[0, column].set_ylim(0.8, 1.35)
    axes[0, 0].set_ylabel("RMSE (mm/h)")
    axes[1, 0].set_ylabel("Bias (mm/h)")
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=4.5,
                      markerfacecolor=SURFACE if style == "hollow" else color,
                      markeredgecolor=color if style == "hollow" else SURFACE, label=METHOD_LABELS[m])
               for m, style, color in PANEL]
    handles.append(Line2D([], [], color=INK, linestyle=(0, (3, 2)), linewidth=0.8, label=METHOD_LABELS["adw"]))
    fig.legend(handles=handles, frameon=False, fontsize=6.0, ncol=len(handles), loc="lower center",
               bbox_to_anchor=(0.5, 0.0))
    save(fig, "lagnbr_fig2_anchors")


# --------------------------------------------------------------------------------------------------
# Figure 3: where the gain comes from, for every anchored method


HIGHLIGHT = {"raw": MUTED, "mgwr": ORDINAL_BLUES[1], "blend_mgwr": ORDINAL_BLUES[3]}


def fig3_strata(comparison: list[dict]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.4), sharex=True)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.84, bottom=0.14, wspace=0.12, hspace=0.12)
    title(fig, "The gain sits in light and moderate rain; heavy rain barely moves",
          "MERGED_OLS_LAGNBR against MERGED_OLS for all 11 methods that read the anchor · thin grey = the other eight")
    x = np.arange(len(STRATA))
    for column, (scheme, heading) in enumerate(SCHEMES):
        top, bottom = axes[0, column], axes[1, column]
        for ax in (top, bottom):
            clean_axis(ax)
            zero_line(ax, vertical=False)
        for method in ANCHORED:
            rows = [one(comparison, scheme=scheme, method=method, stratum=s) for s, _ in STRATA]
            gains = [pct(r["relative_improvement"]) for r in rows]
            bias = [abs(r["Bias_treatment"]) - abs(r["Bias_baseline"]) for r in rows]
            if method in HIGHLIGHT:
                continue
            top.plot(x, gains, color=MUTED, linewidth=0.6, alpha=0.55, zorder=1)
            bottom.plot(x, bias, color=MUTED, linewidth=0.6, alpha=0.55, zorder=1)
        for offset, (method, color) in zip((-0.08, 0.0, 0.08), HIGHLIGHT.items()):
            rows = [one(comparison, scheme=scheme, method=method, stratum=s) for s, _ in STRATA]
            hollow = method == "raw"
            top.plot(x + offset, [pct(r["relative_improvement"]) for r in rows], color=color, linewidth=1.0,
                     linestyle=(0, (2, 1.5)) if hollow else "-", zorder=2)
            bottom.plot(x + offset, [abs(r["Bias_treatment"]) - abs(r["Bias_baseline"]) for r in rows], color=color,
                        linewidth=1.0, linestyle=(0, (2, 1.5)) if hollow else "-", marker="o", markersize=3.2,
                        markerfacecolor=SURFACE if hollow else color, markeredgecolor=color if hollow else SURFACE,
                        zorder=3)
            for xi, r in zip(x + offset, rows):
                ci_dot(top, xi, pct(r["relative_improvement"]), pct(r["ci_low"]), pct(r["ci_high"]), color,
                       vertical=True, hollow=hollow, size=3.6)
        top.set_title(heading, loc="left")
        bottom.set_xticks(x, [label for _, label in STRATA], fontsize=5.8)
    axes[0, 0].set_ylabel("RMSE improvement (%)")
    axes[1, 0].set_ylabel("Change in |bias| (mm/h)\nbelow zero = less biased")
    handles = [Line2D([], [], color=color, linewidth=1.0, linestyle=(0, (2, 1.5)) if m == "raw" else "-",
                      marker="o", markersize=4, markerfacecolor=SURFACE if m == "raw" else color,
                      markeredgecolor=color if m == "raw" else SURFACE, label=METHOD_LABELS[m])
               for m, color in HIGHLIGHT.items()]
    handles.append(Line2D([], [], color=MUTED, linewidth=0.6, alpha=0.7, label="Other eight methods"))
    fig.legend(handles=handles, frameon=False, fontsize=6.2, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 0.0))
    save(fig, "lagnbr_fig3_gain_by_intensity")


# --------------------------------------------------------------------------------------------------
# Figure 4: both merged anchors against gauge-only ADW, with detection skill


DETECTION_METHODS = [("mgwr", "MGWR"), ("blend_mgwr", "Blend")]
THRESHOLDS = ["0.1", "2.5", "8.0", "16.0"]


def fig4_vs_adw(vs_adw: list[dict], detection: list[dict]) -> None:
    fig = plt.figure(figsize=(7.0, 4.7))
    grid = fig.add_gridspec(2, len(STRATA) + 1, width_ratios=[1] * len(STRATA) + [1.25], left=0.155, right=0.915,
                            top=0.82, bottom=0.14, wspace=0.14, hspace=0.36)
    title(fig, "Against gauge-only ADW the new anchor wins more often, most clearly in moderate rain",
          "RMSE improvement over ADW (%) with 95% day-block interval · right: CSI minus ADW's CSI at each threshold "
          "(mm/h); Blend = MGWR + ADW blend")
    positions = np.arange(len(VS_ADW_METHODS))[::-1]
    for row_index, (scheme, heading) in enumerate(SCHEMES):
        for column, (stratum, label) in enumerate(STRATA):
            ax = fig.add_subplot(grid[row_index, column])
            clean_axis(ax, "x")
            zero_line(ax)
            for shift, product, color in ((0.16, OLD, OLD_COLOR), (-0.16, NEW, NEW_COLOR)):
                for y, method in zip(positions, VS_ADW_METHODS):
                    r = one(vs_adw, scheme=scheme, product=product, method=method, stratum=stratum)
                    ci_dot(ax, pct(r["relative_improvement"]), y + shift, pct(r["ci_low"]), pct(r["ci_high"]),
                           color, size=3.2)
            ax.set_ylim(-0.6, len(VS_ADW_METHODS) - 0.4)
            ax.set_yticks(positions, [METHOD_LABELS[m] for m in VS_ADW_METHODS] if column == 0 else [], fontsize=5.4)
            ax.tick_params(labelsize=5.0)
            ax.set_title(label.replace("\n", " "), loc="left", fontsize=6.3)
            if column == 0:
                ax.text(-0.95, 1.13, heading, transform=ax.transAxes, fontsize=7.2, fontweight="semibold", va="bottom")
            if row_index == 1 and column == 2:
                ax.set_xlabel("RMSE improvement over ADW (%)")

        ax = fig.add_subplot(grid[row_index, len(STRATA)])
        clean_axis(ax, "x")
        zero_line(ax)
        labels, y = [], 0
        for threshold in THRESHOLDS[::-1]:
            for method, short in DETECTION_METHODS[::-1]:
                for shift, product, color in ((0.16, OLD, OLD_COLOR), (-0.16, NEW, NEW_COLOR)):
                    r = one(detection, scheme=scheme, product=product, method=method, threshold=float(threshold),
                            metric="CSI")
                    ci_dot(ax, r["delta"], y + shift, r["ci_low"], r["ci_high"], color, size=3.0)
                labels.append(f"{short} · {float(threshold):g}")
                y += 1
        ax.set_yticks(range(len(labels)), labels, fontsize=4.9)
        ax.yaxis.tick_right()
        ax.tick_params(labelsize=5.0)
        ax.set_ylim(-0.6, len(labels) - 0.4)
        ax.set_title("Detection: ΔCSI vs ADW", loc="left", fontsize=6.3)
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=color,
                      markeredgecolor=SURFACE, label=label)
               for color, label in ((OLD_COLOR, "On MERGED_OLS"), (NEW_COLOR, "On MERGED_OLS_LAGNBR"))]
    fig.legend(handles=handles, frameon=False, fontsize=6.4, ncol=2, loc="lower center", bbox_to_anchor=(0.55, 0.0))
    save(fig, "lagnbr_fig4_vs_adw")


# --------------------------------------------------------------------------------------------------
# Figure 5: what the banded fit learned


def band_label(group: int) -> str:
    if group == 0:
        return "No product wet"
    agreement, envelope = divmod(group - 1, 3)
    return f"{agreement + 1} wet · max {ENVELOPE[envelope]}"


TERM_GROUPS = [("", "At hour t", SERIES[0]), ("_lag1", "t − 1 h", SERIES[1]), ("_lead1", "t + 1 h", SERIES[2]),
               ("_nbr8", "Nearest 8 stations", SERIES[3])]


def fig5_coefficients(coefficients: list[dict]) -> None:
    rows = [r for r in coefficients if r["scheme"] == "balanced_spatial"]
    get = lambda g, t, c: one(rows, group=float(g), term=t)[c]
    matrix = np.array([[get(g, t, "coefficient_mean") for t in TERMS] for g in range(10)])
    spread = np.array([[get(g, t, "coefficient_max") - get(g, t, "coefficient_min") for t in TERMS] for g in range(10)])
    fallback = [get(g, "intercept", "folds_fell_back") > 0 for g in range(10)]
    cells = np.array([get(g, "intercept", "n_cell_mean") for g in range(10)])

    fig = plt.figure(figsize=(7.0, 4.1))
    grid = fig.add_gridspec(1, 3, width_ratios=[13, 2.2, 3.4], left=0.175, right=0.99, top=0.8, bottom=0.15,
                            wspace=0.06)
    title(fig, "What the banded fit learned: GPM's neighbourhood mean carries the most weight in every band",
          "Fold-mean coefficient (small: ± half the fold-to-fold range) · † = band took the pooled fit in every fold")
    ax = fig.add_subplot(grid[0, 0])
    limit = np.nanmax(np.abs(matrix[:, 1:]))
    ax.imshow(matrix, cmap=DIVERGING, norm=TwoSlopeNorm(0, -limit, limit), aspect="auto")
    for (i, j), v in np.ndenumerate(matrix):
        shade = abs(v) / limit > 0.55
        ax.text(j, i - 0.13, f"{v:.2f}", ha="center", va="center", fontsize=4.9,
                color=SURFACE if shade else INK_SECONDARY)
        ax.text(j, i + 0.24, f"±{spread[i, j] / 2:.2f}", ha="center", va="center", fontsize=3.7,
                color=SURFACE if shade else MUTED)
    ax.set_xticks(range(len(TERMS)), ["Int.", "FY4B", "GPM", "GSMaP", "FY4B", "GPM", "GSMaP", "FY4B", "GPM",
                                      "GSMaP", "FY4B", "GPM", "GSMaP"], fontsize=5.2)
    for start, label in ((1, "at hour t"), (4, "t − 1 h"), (7, "t + 1 h"), (10, "nearest 8 stations")):
        ax.text(start + 1, -0.75, label, ha="center", va="bottom", fontsize=5.8, color=INK_SECONDARY)
        ax.axvline(start - 0.5, color=SURFACE, linewidth=1.6)
    ax.set_yticks(range(10), [f"{band_label(g)}{' †' if fallback[g] else ''}" for g in range(10)], fontsize=5.6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    bars = fig.add_subplot(grid[0, 1])
    bars.barh(range(10), cells, height=0.6, color=[MUTED if f else NEW_COLOR for f in fallback], zorder=2)
    bars.set_xscale("log")
    bars.set_ylim(9.5, -0.5)
    bars.set_yticks([])
    clean_axis(bars, "x")
    bars.tick_params(labelsize=4.8)
    bars.set_title("Training cells\nper fold", loc="left", fontsize=5.8, color=INK_SECONDARY)
    for g, n in enumerate(cells):
        bars.text(n * 1.15, g, f"{n / 1e3:,.0f}k", va="center", fontsize=4.4, color=INK_SECONDARY)
    bars.set_xlim(500, 2e7)

    sums = fig.add_subplot(grid[0, 2])
    clean_axis(sums, "x")
    zero_line(sums)
    for offset, (suffix, label, color) in zip((-0.24, -0.08, 0.08, 0.24), TERM_GROUPS):
        values = [sum(get(g, f"{p}{suffix}", "coefficient_mean") for p in ("fy4b", "gpm", "gsmap")) for g in range(10)]
        sums.plot(values, np.arange(10) + offset, linestyle="none", marker="o", markersize=3.3, markerfacecolor=color,
                  markeredgecolor=SURFACE, markeredgewidth=0.4, label=label, zorder=3)
    sums.set_ylim(9.5, -0.5)
    sums.set_xlim(-0.12, 0.58)
    sums.set_yticks([])
    sums.tick_params(labelsize=4.8)
    sums.set_title("Summed weight\nper term group", loc="left", fontsize=5.8, color=INK_SECONDARY)
    fig.legend(*sums.get_legend_handles_labels(), frameon=False, fontsize=6.0, ncol=4, loc="lower center",
               bbox_to_anchor=(0.6, 0.0))
    save(fig, "lagnbr_fig5_coefficients")


# --------------------------------------------------------------------------------------------------
# Report tables: rounded strings, so the Typst file only lays them out.


def fmt(x, digits=3, signed=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:+.{digits}f}" if signed else f"{x:.{digits}f}"


def write_table(name: str, header: list[str], rows: list[list[str]]) -> None:
    """One report table as JSON records of pre-formatted strings (JSON because the repository ignores *.csv)."""
    path = REPORT_DIR / "tables" / f"{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    records = [dict(zip(header, row)) for row in rows]
    path.write_text(json.dumps(records, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"  wrote tables/{name}.json")


def paired_cells(row) -> list[str]:
    return [fmt(row["RMSE_baseline"], 4), fmt(row["RMSE_treatment"], 4),
            fmt(pct(row["relative_improvement"]), 2, True), fmt(pct(row["ci_low"]), 2, True),
            fmt(pct(row["ci_high"]), 2, True), fmt(row["pvalue_holm"], 3)]


def report_tables(summary, comparison, anchors, anchor_tests, vs_adw, detection, coefficients, claim, coverage,
                  invariance, shared, run_dir) -> None:
    header = ["scheme", "method", "stratum", "n", "RMSE_ols", "RMSE_new", "gain", "ci_low", "ci_high", "p_holm"]
    write_table("lagnbr_vs_ols", header, [
        [r["scheme"], r["method"], r["stratum"], f"{int(r['n']):,}"] + paired_cells(r) for r in comparison])

    rows = []
    for scheme, _ in SCHEMES:
        for anchor in ANCHORS:
            for method in ("raw", "residual_gwr", "mgwr", "auto", "blend_mgwr", "blend_agrenv_mgwr"):
                test = one(anchor_tests, scheme=scheme, product=anchor, method=method)
                rmse = {s: one(anchors, scheme=scheme, product=anchor, method=method, stratum=s)["RMSE"]
                        for s in ("all", "moderate", "heavy")}
                rows.append([scheme, anchor, method, fmt(rmse["all"], 4), fmt(rmse["moderate"], 3),
                             fmt(rmse["heavy"], 2),
                             "–" if test is None else fmt(pct(test["relative_improvement"]), 2, True),
                             "–" if test is None else fmt(pct(test["ci_low"]), 2, True),
                             "–" if test is None else fmt(pct(test["ci_high"]), 2, True)])
        adw = {s: one(anchors, scheme=scheme, product="gauge_only", method="adw", stratum=s)["RMSE"]
               for s in ("all", "moderate", "heavy")}
        rows.append([scheme, "gauge_only", "adw", fmt(adw["all"], 4), fmt(adw["moderate"], 3), fmt(adw["heavy"], 2),
                     "–", "–", "–"])
    write_table("lagnbr_anchors", ["scheme", "anchor", "method", "RMSE", "RMSE_moderate", "RMSE_heavy", "vs_ols",
                                   "ci_low", "ci_high"], rows)

    write_table("lagnbr_vs_adw", ["scheme", "product", "method", "stratum", "RMSE_adw", "RMSE_method", "gain",
                                  "ci_low", "ci_high", "p_holm"],
                [[r["scheme"], r["product"], r["method"], r["stratum"]] + paired_cells(r) for r in vs_adw])

    write_table("lagnbr_detection", ["scheme", "product", "method", "threshold", "metric", "adw", "method_value",
                                     "delta", "ci_low", "ci_high"],
                [[r["scheme"], r["product"], r["method"], fmt(r["threshold"], 1), r["metric"], fmt(r["value_adw"], 3),
                  fmt(r["value_method"], 3), fmt(r["delta"], 3, True), fmt(r["ci_low"], 3, True),
                  fmt(r["ci_high"], 3, True)] for r in detection])

    rows = []
    for scheme, _ in SCHEMES:
        for product in ANCHORS:
            for method in ALL_METHODS:
                r = one(summary, scheme=scheme, product=product, method=method)
                if r is not None:
                    rows.append([scheme, product, method, f"{int(r['n']):,}", fmt(r["coverage"], 3), fmt(r["RMSE"], 4),
                                 fmt(r["Bias"], 4, True), fmt(r["r"], 3),
                                 fmt(pct(r["RMSE_improvement_vs_raw"]), 1, True),
                                 fmt(pct(r["RMSE_improvement_vs_adw"]), 1, True)])
    write_table("lagnbr_method_summary", ["scheme", "product", "method", "n", "coverage", "RMSE", "Bias", "r", "vs_raw",
                                          "vs_adw"], rows)

    rows = []
    for g in range(10):
        members = [r for r in coefficients if r["scheme"] == "balanced_spatial" and r["group"] == float(g)]
        by_term = {r["term"]: r for r in members}
        group_sum = lambda suffix: sum(by_term[f"{p}{suffix}"]["coefficient_mean"] for p in ("fy4b", "gpm", "gsmap"))
        first = by_term["intercept"]
        rows.append([str(g), band_label(g), f"{first['n_cell_mean']:,.0f}", f"{int(first['folds_fell_back'])}/5",
                     fmt(first["coefficient_mean"], 3), fmt(group_sum(""), 3), fmt(group_sum("_lag1"), 3),
                     fmt(group_sum("_lead1"), 3), fmt(group_sum("_nbr8"), 3)])
    write_table("lagnbr_bands", ["band", "label", "cells", "fell_back", "intercept", "sum_t", "sum_lag", "sum_lead",
                                 "sum_nbr"], rows)

    write_table("lagnbr_claim", ["product", "RMSE_residual_gwr", "overall_gain", "paired_significant", "heavy_gain",
                                 "heavy_ci_low", "moderate_degradation", "coverage", "coverage_ok", "events_ok",
                                 "supported"],
                [[r["product"], fmt(r["RMSE_residual_gwr"], 4), fmt(pct(r["overall_relative_improvement"]), 2, True),
                  "yes" if r["paired_significant"] else "no", fmt(pct(r["heavy_relative_improvement"]), 2, True),
                  fmt(pct(r["heavy_relative_ci_low"]), 2, True), fmt(pct(r["moderate_relative_degradation"]), 2, True),
                  fmt(r["common_coverage"], 3), "yes" if r["coverage_acceptable"] else "no",
                  "yes" if r["event_not_degraded"] else "no", "yes" if r["product_supported"] else "no"]
                 for r in claim])

    write_table("lagnbr_coverage", ["scheme", "product", "fold", "status", "coverage", "covariates"],
                [[r["scheme"], r["product"], str(int(r["fold"])), r["status"], fmt(r["prediction_coverage"], 3),
                  r["covariate_variables"]] for r in coverage
                 if r["method"] == "mgwr" and r["scheme"] == "balanced_spatial"])

    write_table("lagnbr_invariance", ["file", "baseline_rows", "identical", "changed"],
                [[r["file"], str(int(r["baseline_rows"])), str(int(r["identical"])), str(int(r["changed"]))]
                 for r in invariance])
    write_table("lagnbr_shared_cells", ["scheme", "comparison", "cells"],
                [[r["scheme"], r["comparison"], f"{int(r['shared_cells']):,}"] for r in shared])

    scope = read_rows(run_dir / "benchmark_scope.csv")
    commit = next((r["value"] for r in scope if r["key"] == "git_commit"), "unknown")
    write_table("lagnbr_provenance", ["key", "value"], [["run", run_dir.name], ["git_commit", str(commit)[:7]]])


def main() -> None:
    run_dir = Path((RUN_DIR / "source_run.txt").read_text(encoding="utf-8").strip())
    summary = read_rows(RUN_DIR / "method_summary.csv")
    comparison = read_rows(RUN_DIR / "lagnbr_vs_ols.csv")
    anchors = read_rows(RUN_DIR / "anchor_metrics.csv")
    anchor_tests = read_rows(RUN_DIR / "anchor_comparison.csv")
    vs_adw = read_rows(RUN_DIR / "vs_adw.csv")
    detection = read_rows(RUN_DIR / "detection_skill.csv")
    coefficients = read_rows(RUN_DIR / "grouped_coefficients_summary.csv")
    claim = read_rows(RUN_DIR / "claim_assessment.csv")
    coverage = read_rows(RUN_DIR / "coverage_note.csv")
    invariance = read_rows(RUN_DIR / "invariance.csv")
    shared = read_rows(RUN_DIR / "shared_cells.csv")

    fig1_gain(comparison)
    fig2_anchors(anchors)
    fig3_strata(comparison)
    fig4_vs_adw(vs_adw, detection)
    fig5_coefficients(coefficients)

    (REPORT_DIR / "figures").mkdir(parents=True, exist_ok=True)
    for path in sorted(FIG_DIR.glob("lagnbr_fig*.png")):
        shutil.copy2(path, REPORT_DIR / "figures" / path.name)
    report_tables(summary, comparison, anchors, anchor_tests, vs_adw, detection, coefficients, claim, coverage,
                  invariance, shared, run_dir)
    print(f"Copied figures and tables into {REPORT_DIR}")


if __name__ == "__main__":
    main()
