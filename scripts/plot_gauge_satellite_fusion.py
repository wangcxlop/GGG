"""Draw the gauge-satellite fusion figures and the report tables.

Reads the CSVs written by `scripts/run_gauge_satellite_fusion.jl` and writes PNG (300 dpi) and PDF figures
to `output/gauge_satellite_fusion/figures/`, then copies the figures and a set of rounded report tables
into `Interim_results/` so the Typst report compiles from that folder alone:

    py -3.13 scripts/plot_gauge_satellite_fusion.py

Colour carries one job per figure. Anchors (FY4B, GPM, GSMaP and the two merges) sit on axes, never on
colour, so the five of them never compete for categorical slots. Where colour marks how much the gauges
and satellites are combined, it is an ordinal ramp: hollow grey for the raw satellite field, a dashed ink
line for gauge-only ADW, then two steps of the blue ramp for the in-fold GWR-family choice (`auto`) and
MGWR with the agreement blend. Improvements use a blue (better) / red (worse) diverging ramp with a grey
centre.
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
RUN_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "output" / "gauge_satellite_fusion"
FIG_DIR = RUN_DIR / "figures"
REPORT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "Interim_results"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
ORDINAL_BLUES = ["#86b6ef", "#5598e7", "#256abf", "#104281"]
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # validated categorical slots 1-4 (adjacent use)
DIVERGING = LinearSegmentedColormap.from_list(
    "better_worse", ["#b52f2f", "#e34948", "#ef8d8c", "#f0efec", "#86b6ef", "#3987e5", "#1c5cab"])

SCHEME = "balanced_spatial"
ANCHORS = ["FY4B", "GPM", "GSMaP", "MERGED_MEAN", "MERGED_OLS"]
ANCHOR_LABELS = {"FY4B": "FY4B", "GPM": "GPM", "GSMaP": "GSMaP", "MERGED_MEAN": "Merged\n(equal)",
                 "MERGED_OLS": "Merged\n(OLS)"}
FUSION_METHODS = ["residual_gwr", "mixed_gwr", "mgwr", "blend_residual_gwr", "blend_mixed_gwr", "blend_mgwr",
                  "blend_agrenv_residual_gwr", "blend_agrenv_mixed_gwr", "blend_agrenv_mgwr"]
METHOD_LABELS = {
    "raw": "Raw satellite", "adw": "ADW (gauges only)", "idw": "IDW", "tps": "TPS", "gwr": "GWR (gauges only)",
    "auto": "GWR family, chosen in-fold", "residual_gwr": "Residual GWR", "mixed_gwr": "Mixed GWR", "mgwr": "MGWR",
    "blend_residual_gwr": "Residual GWR + ADW blend", "blend_mixed_gwr": "Mixed GWR + ADW blend",
    "blend_mgwr": "MGWR + ADW blend", "blend_agrenv_residual_gwr": "Residual GWR + agreement blend",
    "blend_agrenv_mixed_gwr": "Mixed GWR + agreement blend", "blend_agrenv_mgwr": "MGWR + agreement blend",
}
# The integration ladder, in order. Raw and ADW are drawn as references, not as coloured series.
LADDER = [("auto", ORDINAL_BLUES[1]), ("blend_agrenv_mgwr", ORDINAL_BLUES[3])]
INTENSITY = [("no_rain", "Dry\n< 0.1"), ("light", "Light\n0.1–2.5"), ("moderate", "Moderate\n2.5–8"),
             ("heavy", "Heavy\n≥ 8")]
DISTANCE = [("0_20", "< 20 km"), ("20_50", "20–50 km"), ("50_100", "50–100 km")]
THRESHOLDS = [0.1, 2.5, 8.0, 16.0]

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
                if key in ("key", "value"):   # provenance: a commit hash must not become a number
                    parsed[key] = value
                    continue
                if value in ("true", "false"):
                    parsed[key] = value == "true"
                    continue
                if key in ("level", "product", "method", "group", "stratum", "baseline", "scheme", "metric", "axis"):
                    parsed[key] = value
                    continue
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


def value(rows, column, **conditions) -> float:
    row = one(rows, **conditions)
    return float("nan") if row is None else float(row[column])


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
    # A fixed physical gap, so short and tall figures space the subtitle the same way.
    fig.text(0.02, top - 0.24 / fig.get_figheight(), subheading, ha="left", va="top", fontsize=6.3,
             color=INK_SECONDARY)


def ladder_legend(fig: plt.Figure, raw: bool = True, **kwargs) -> None:
    handles = []
    if raw:
        handles.append(Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=SURFACE,
                              markeredgecolor=MUTED, label=METHOD_LABELS["raw"]))
    handles.append(Line2D([], [], color=INK, linestyle=(0, (3, 2)), linewidth=0.9, label=METHOD_LABELS["adw"]))
    handles += [Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=color,
                       markeredgecolor=SURFACE, label=METHOD_LABELS[method]) for method, color in LADDER]
    fig.legend(handles=handles, frameon=False, fontsize=6.6, ncol=len(handles), **kwargs)


def pct(x: float) -> float:
    return 100.0 * x


# --------------------------------------------------------------------------------------------------
# Figure 1: the benefit ladder


BENEFIT_PANELS = [("RMSE", "RMSE (mm/h) · lower is better"), ("r", "Correlation r · higher is better"),
                  ("CSI_2.5", "CSI at 2.5 mm/h · higher is better"), ("CSI_8.0", "CSI at 8 mm/h · higher is better")]


def fig1_benefit(summary: list[dict]) -> None:
    fig, axes = plt.subplots(1, len(BENEFIT_PANELS), figsize=(10.0, 3.3))
    fig.subplots_adjust(left=0.05, right=0.99, top=0.76, bottom=0.2, wspace=0.28)
    x = np.arange(len(ANCHORS))
    offsets = {"raw": -0.22, "auto": 0.0, "blend_agrenv_mgwr": 0.22}
    for ax, (column, label) in zip(axes, BENEFIT_PANELS):
        clean_axis(ax)
        for index, anchor in enumerate(ANCHORS):
            reference = value(summary, column, scheme=SCHEME, product=anchor, method="adw")
            ax.hlines(reference, index - 0.36, index + 0.36, color=INK, linestyle=(0, (3, 2)), linewidth=0.9, zorder=2)
        raw = [value(summary, column, scheme=SCHEME, product=a, method="raw") for a in ANCHORS]
        ax.scatter(x + offsets["raw"], raw, s=24, facecolors=SURFACE, edgecolors=MUTED, linewidths=1.0, zorder=3)
        for method, color in LADDER:
            values = [value(summary, column, scheme=SCHEME, product=a, method=method) for a in ANCHORS]
            ax.scatter(x + offsets[method], values, s=26, color=color, edgecolors=SURFACE, linewidths=0.8, zorder=4)
        ax.set_xticks(x)
        ax.set_xticklabels([ANCHOR_LABELS[a] for a in ANCHORS], fontsize=6)
        ax.set_title(label, loc="left", pad=3)
    title(fig, "Combining gauges with satellite rainfall: skill at gauges held out of the fit",
          "Balanced spatial 5-fold cross-validation, 237 gauges, 13,471 hours (2022–2024). Each anchor is scored on "
          "its own evaluation mask; the dashed line is gauge-only ADW on that mask.", top=0.97)
    ladder_legend(fig, loc="lower center", bbox_to_anchor=(0.5, 0.0))
    save(fig, "fusion_fig1_benefit")


# --------------------------------------------------------------------------------------------------
# Figure 2: interpolation methods


def fig2_interpolators(ranking: list[dict], summary: list[dict]) -> None:
    fig = plt.figure(figsize=(10.0, 3.2))
    strata = [("overall", "All hours"), ("moderate", "Moderate hours (2.5–8 mm/h)"), ("heavy", "Heavy hours (≥ 8 mm/h)")]
    pairs = [("idw", "adw"), ("tps", "adw"), ("gwr", "adw")]
    width = 0.19
    for panel, (stratum, label) in enumerate(strata):
        ax = fig.add_axes([0.13 + panel * 0.185, 0.2, 0.155, 0.55])
        clean_axis(ax, "x")
        ax.axvline(0, color=AXIS, linewidth=0.8, zorder=1)
        for row_index, (other, reference) in enumerate(pairs):
            # The ranking tables store `method - baseline`; express every row as "other vs ADW".
            if other == "gwr":
                row = one(ranking, method="gwr", baseline="adw", stratum=stratum)
                sign = 1.0
            else:
                row = one(ranking, method="adw", baseline=other, stratum=stratum)
                sign = -1.0
            if row is None:
                continue
            base = row["RMSE_gwr"] if sign < 0 else row["RMSE_baseline"]
            centre = pct(sign * row["delta_RMSE"] / base)
            lo, hi = sorted([pct(sign * row["ci_low"] / base), pct(sign * row["ci_high"] / base)])
            y = len(pairs) - 1 - row_index
            ax.hlines(y, lo, hi, color=INK_SECONDARY, linewidth=1.2, zorder=2)
            ax.scatter([centre], [y], s=24, color=INK, edgecolors=SURFACE, linewidths=0.8, zorder=3)
        ax.set_yticks(range(len(pairs)))
        ax.set_yticklabels([METHOD_LABELS[p[0]] for p in reversed(pairs)] if panel == 0 else [], fontsize=6)
        ax.set_title(label, loc="left", pad=3)
        ax.set_xlabel("RMSE change vs ADW (%) · right is better", fontsize=6)
    # Grouped dots rather than lines: IDW and ADW agree to the third decimal, and two lines drawn on
    # top of each other would hide one of them.
    ax = fig.add_axes([0.71, 0.2, 0.27, 0.55])
    clean_axis(ax)
    interpolators = ["idw", "adw", "tps", "gwr"]
    for index, method in enumerate(interpolators):
        values = [value(summary, f"CSI_{t}", scheme=SCHEME, product="GPM", method=method) for t in THRESHOLDS]
        ax.scatter(np.arange(len(THRESHOLDS)) + (index - 1.5) * 0.17, values, s=20, color=SERIES[index],
                   edgecolors=SURFACE, linewidths=0.8, zorder=3)
    ax.set_xticks(range(len(THRESHOLDS)))
    ax.set_xticklabels([f"≥ {t:g}" for t in THRESHOLDS])
    ax.set_xlabel("Event threshold (mm/h)", fontsize=6)
    ax.set_title("Critical success index (higher is better)", loc="left", pad=3)
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=SERIES[i],
                      markeredgecolor=SURFACE, label=METHOD_LABELS[m]) for i, m in enumerate(interpolators)]
    fig.legend(handles=handles, frameon=False, fontsize=6.4, ncol=4, loc="lower right", bbox_to_anchor=(0.99, 0.0))
    title(fig, "Gauge-only interpolation: ADW and IDW lead, TPS and GWR trail",
          "Held-out gauges, balanced spatial CV. Whiskers: 95% day-block bootstrap interval of the RMSE difference "
          "against ADW. Interpolators ignore the satellite, so GPM's evaluation mask is used.", top=0.97)
    save(fig, "fusion_fig2_interpolators")


# --------------------------------------------------------------------------------------------------
# Figure 3: fusion method x anchor


def fig3_heatmap(fusion: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0))
    fig.subplots_adjust(left=0.17, right=0.93, top=0.8, bottom=0.12, wspace=0.08)
    norm = TwoSlopeNorm(vmin=-16, vcenter=0, vmax=16)
    for ax, (stratum, label) in zip(axes, [("overall", "All hours"), ("heavy", "Heavy hours (≥ 8 mm/h)")]):
        grid = np.full((len(FUSION_METHODS), len(ANCHORS)), np.nan)
        for i, method in enumerate(FUSION_METHODS):
            for j, anchor in enumerate(ANCHORS):
                row = one(fusion, stratum=stratum, product=anchor, method=method)
                if row is None:
                    continue
                grid[i, j] = pct(row["relative_improvement"])
                significant = row["pvalue_holm_family"] < 0.05 and np.sign(row["relative_ci_low"]) == np.sign(
                    row["relative_ci_high"])
                text = f"{grid[i, j]:+.1f}" + ("*" if significant else "")
                shade = abs(grid[i, j]) > 9
                ax.text(j, i, text, ha="center", va="center", fontsize=6,
                        color=SURFACE if shade else INK, fontweight="bold" if significant else "normal")
        image = ax.imshow(np.clip(grid, -16, 16), cmap=DIVERGING, norm=norm, aspect="auto")
        ax.set_xticks(range(len(ANCHORS)))
        ax.set_xticklabels([ANCHOR_LABELS[a] for a in ANCHORS], fontsize=6)
        ax.set_yticks(range(len(FUSION_METHODS)))
        ax.set_yticklabels([METHOD_LABELS[m] for m in FUSION_METHODS] if ax is axes[0] else [], fontsize=6)
        ax.tick_params(length=0, pad=2)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks(np.arange(-0.5, len(ANCHORS)), minor=True)
        ax.set_yticks(np.arange(-0.5, len(FUSION_METHODS)), minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=2)
        ax.tick_params(which="minor", length=0)
        ax.set_title(label, loc="left", pad=4)
    bar = fig.colorbar(image, ax=axes, fraction=0.02, pad=0.015)
    bar.set_label("RMSE improvement vs ADW (%)", fontsize=6)
    bar.ax.tick_params(labelsize=6, length=2)
    bar.outline.set_visible(False)
    title(fig, "Which fusion method, on which anchor, beats gauge-only ADW",
          "Relative RMSE improvement over ADW at held-out gauges (blue better, red worse; colour clipped at ±16%). "
          "* and bold: 95% bootstrap interval excludes 0 and Holm-adjusted p < 0.05 across every cell of the panel.",
          top=0.97)
    save(fig, "fusion_fig3_methods_by_anchor")


# --------------------------------------------------------------------------------------------------
# Figure 4: where the satellite pays


def fig4_where(stratified: list[dict], detection: list[dict]) -> None:
    fig, axes = plt.subplots(3, len(ANCHORS), figsize=(10.0, 7.0), sharey="row")
    fig.subplots_adjust(left=0.07, right=0.99, top=0.84, bottom=0.09, wspace=0.08, hspace=0.72)
    offsets = {"raw": -0.24, "auto": 0.0, "blend_agrenv_mgwr": 0.24}
    for column, anchor in enumerate(ANCHORS):
        for row_index, (group, levels) in enumerate([("rain_intensity", INTENSITY), ("nearest_train_km", DISTANCE)]):
            ax = axes[row_index, column]
            clean_axis(ax)
            ax.axhline(0, color=INK, linestyle=(0, (3, 2)), linewidth=0.9, zorder=1)
            x = np.arange(len(levels))
            for method, color in LADDER:
                values = [pct(value(stratified, "RMSE_improvement_vs_adw", product=anchor, method=method, group=group,
                                    level=level)) for level, _ in levels]
                ax.scatter(x + offsets[method], values, s=22, color=color, edgecolors=SURFACE, linewidths=0.8, zorder=3)
            ax.set_xticks(x)
            ax.set_xticklabels([label for _, label in levels], fontsize=5.6)
            if column == 0:
                ax.set_ylabel("RMSE change vs ADW (%)", fontsize=6)
            if row_index == 0:
                ax.set_title(ANCHOR_LABELS[anchor].replace("\n", " "), loc="left", pad=3)
        ax = axes[2, column]
        clean_axis(ax)
        ax.axhline(0, color=INK, linestyle=(0, (3, 2)), linewidth=0.9, zorder=1)
        x = np.arange(len(THRESHOLDS))
        for method, color in [("raw", None)] + LADDER:
            rows = [one(detection, product=anchor, method=method, threshold=t, metric="CSI") for t in THRESHOLDS]
            centre = [row["delta"] if row else np.nan for row in rows]
            lo = [row["ci_low"] if row else np.nan for row in rows]
            hi = [row["ci_high"] if row else np.nan for row in rows]
            if method == "raw":
                ax.vlines(x + offsets[method], lo, hi, color=MUTED, linewidth=1.0, zorder=2)
                ax.scatter(x + offsets[method], centre, s=20, facecolors=SURFACE, edgecolors=MUTED, linewidths=1.0,
                           zorder=3)
            else:
                ax.vlines(x + offsets[method], lo, hi, color=color, linewidth=1.0, zorder=2)
                ax.scatter(x + offsets[method], centre, s=22, color=color, edgecolors=SURFACE, linewidths=0.8, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([f"≥ {t:g}" for t in THRESHOLDS], fontsize=5.8)
        ax.set_xlabel("Event threshold (mm/h)", fontsize=6)
        if column == 0:
            ax.set_ylabel("CSI change vs ADW", fontsize=6)
    sections = ["RMSE change vs ADW by gauge rain intensity (mm/h)",
                "RMSE change vs ADW by distance from the held-out gauge to the nearest training gauge",
                "CSI change vs ADW, with 95% day-block bootstrap intervals"]
    for row_index, text in enumerate(sections):
        top = axes[row_index, 0].get_position().y1
        # Row 0 carries the anchor names as panel titles, so its section label sits above them.
        fig.text(0.07, top + (0.045 if row_index == 0 else 0.018), text, fontsize=7.2, color=INK,
                 fontweight="semibold")
    fig.suptitle("Where the satellite adds to the gauges: heavy rain, detection, and distance from the network",
                 x=0.02, ha="left", fontsize=9, y=0.985)
    ladder_legend(fig, loc="upper right", bbox_to_anchor=(0.99, 0.955))
    save(fig, "fusion_fig4_where_it_helps")


# --------------------------------------------------------------------------------------------------
# Figure 5: which product


def fig5_products(anchor_metrics: list[dict], comparison: list[dict], coefficients: list[dict]) -> None:
    fig = plt.figure(figsize=(10.0, 6.0))
    panels = [("RMSE", "all", "overall", "RMSE (mm/h)"), ("r", "all", "overall", "Correlation r"),
              ("CSI", "0.1", "event_threshold", "CSI ≥ 0.1 mm/h"), ("CSI", "8.0", "event_threshold", "CSI ≥ 8 mm/h")]
    x = np.arange(len(ANCHORS))
    for index, (column, level, group, label) in enumerate(panels):
        ax = fig.add_axes([0.05 + index * 0.24, 0.6, 0.2, 0.26])
        clean_axis(ax)
        raw = [value(anchor_metrics, column, product=a, method="raw", group=group, level=level) for a in ANCHORS]
        fused = [value(anchor_metrics, column, product=a, method="blend_agrenv_mgwr", group=group, level=level)
                 for a in ANCHORS]
        reference = value(anchor_metrics, column, product="gauge_only", method="adw", group=group, level=level)
        ax.axhline(reference, color=INK, linestyle=(0, (3, 2)), linewidth=0.9, zorder=1)
        ax.bar(x - 0.19, raw, width=0.36, color=AXIS, edgecolor=SURFACE, linewidth=0.8, zorder=2)
        ax.bar(x + 0.19, fused, width=0.36, color=ORDINAL_BLUES[3], edgecolor=SURFACE, linewidth=0.8, zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels([ANCHOR_LABELS[a] for a in ANCHORS], fontsize=5.6)
        ax.set_title(label, loc="left", pad=3)
    fig.text(0.05, 0.9, "Every anchor on the same cells: raw field (grey), after MGWR + agreement blend (blue), "
             "gauge-only ADW (dashed)", fontsize=7)

    ax = fig.add_axes([0.1, 0.16, 0.47, 0.3])
    clean_axis(ax, "x")
    ax.axvline(0, color=INK, linestyle=(0, (3, 2)), linewidth=0.9, zorder=1)
    others = ["FY4B", "GSMaP", "MERGED_MEAN", "MERGED_OLS"]
    methods = [("raw", MUTED, True), ("auto", ORDINAL_BLUES[1], False), ("blend_agrenv_mgwr", ORDINAL_BLUES[3], False)]
    for row_index, anchor in enumerate(others):
        for offset, (method, color, hollow) in zip((-0.22, 0.0, 0.22), methods):
            row = one(comparison, product=anchor, method=method)
            if row is None:
                continue
            y = len(others) - 1 - row_index + offset
            ax.hlines(y, pct(row["relative_ci_low"]), pct(row["relative_ci_high"]), color=color, linewidth=1.1, zorder=2)
            ax.scatter([pct(row["relative_improvement"])], [y], s=22, zorder=3,
                       facecolors=SURFACE if hollow else color, edgecolors=color if hollow else SURFACE,
                       linewidths=1.0 if hollow else 0.8)
    ax.set_yticks(range(len(others)))
    ax.set_yticklabels([ANCHOR_LABELS[a].replace("\n", " ") for a in reversed(others)], fontsize=6)
    ax.set_xlabel("RMSE change vs the same method on GPM (%) · right is better", fontsize=6)
    ax.set_title("Each anchor against GPM, same method on both sides, 95% bootstrap intervals", loc="left", pad=3)
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=SURFACE if h else c,
                      markeredgecolor=c if h else SURFACE, label=METHOD_LABELS[m]) for m, c, h in methods]
    ax.legend(handles=handles, frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(0.0, -0.13), ncol=3)

    ax = fig.add_axes([0.66, 0.16, 0.32, 0.3])
    clean_axis(ax)
    ols = one(coefficients, scheme=SCHEME, product="MERGED_OLS")
    names = [("intercept", "Intercept"), ("beta_fy4b", "FY4B"), ("beta_gpm", "GPM"), ("beta_gsmap", "GSMaP")]
    if ols is not None:
        for index, (key, _) in enumerate(names):
            ax.vlines(index, ols[f"{key}_min"], ols[f"{key}_max"], color=ORDINAL_BLUES[3], linewidth=6, alpha=0.35,
                      zorder=2)
            ax.scatter([index], [ols[f"{key}_mean"]], s=26, color=ORDINAL_BLUES[3], edgecolors=SURFACE, zorder=3)
            ax.annotate(f"{ols[f'{key}_mean']:.2f}", (index, ols[f"{key}_mean"]), xytext=(6, 0),
                        textcoords="offset points", fontsize=6, color=INK_SECONDARY, va="center")
    ax.axhline(0, color=AXIS, linewidth=0.8, zorder=1)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([label for _, label in names], fontsize=6)
    ax.set_title("OLS merge weights: mean (dot) and range over the 5 folds (band)", loc="left", pad=3)
    ax.set_ylabel("Weight (intercept in mm/h)", fontsize=6)
    fig.suptitle("Which satellite product, or should all three be merged?", x=0.02, ha="left", fontsize=9, y=0.985)
    save(fig, "fusion_fig5_products")


# --------------------------------------------------------------------------------------------------
# Figure 6: robustness


def fig6_robustness(summary: list[dict], stratified: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.2), gridspec_kw={"width_ratios": [1, 1, 1.5]})
    fig.subplots_adjust(left=0.05, right=0.99, top=0.76, bottom=0.2, wspace=0.28)
    x = np.arange(len(ANCHORS))
    for ax, (scheme, label) in zip(axes[:2], [("balanced_spatial", "Spatial CV (held-out regions)"),
                                              ("random", "Random CV (held-out gauges among neighbours)")]):
        clean_axis(ax)
        for index, anchor in enumerate(ANCHORS):
            reference = value(summary, "RMSE", scheme=scheme, product=anchor, method="adw")
            ax.hlines(reference, index - 0.36, index + 0.36, color=INK, linestyle=(0, (3, 2)), linewidth=0.9, zorder=2)
        for offset, (method, color) in zip((-0.12, 0.12), LADDER):
            values = [value(summary, "RMSE", scheme=scheme, product=a, method=method) for a in ANCHORS]
            ax.scatter(x + offset, values, s=24, color=color, edgecolors=SURFACE, linewidths=0.8, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([ANCHOR_LABELS[a] for a in ANCHORS], fontsize=5.8)
        ax.set_title(f"RMSE · {label}", loc="left", pad=3)
    ax = axes[2]
    clean_axis(ax)
    ax.axhline(0, color=INK, linestyle=(0, (3, 2)), linewidth=0.9, zorder=1)
    # The window's last 08-08 day ends at 08:00 on 1 January 2025, so a handful of hours carry the label
    # 2025; they are not a year and are left out.
    years = sorted({row["level"] for row in stratified if row["group"] == "year" and row["level"] != "2025"})
    width = 0.8 / len(years)
    for year_index, year in enumerate(years):
        shade = ORDINAL_BLUES[1 + year_index % (len(ORDINAL_BLUES) - 1)]
        values = [pct(value(stratified, "RMSE_improvement_vs_adw", product=a, method="blend_agrenv_mgwr", group="year",
                            level=year)) for a in ANCHORS]
        ax.bar(x - 0.4 + width * (year_index + 0.5), values, width=width * 0.9, color=shade, edgecolor=SURFACE,
               linewidth=0.6, zorder=2, label=year)
    ax.set_xticks(x)
    ax.set_xticklabels([ANCHOR_LABELS[a] for a in ANCHORS], fontsize=5.8)
    ax.set_title("MGWR + agreement blend vs ADW by year (RMSE change, %)", loc="left", pad=3)
    ax.legend(frameon=False, fontsize=6, ncol=len(years), loc="lower right", bbox_to_anchor=(1.0, -0.3))
    title(fig, "Robustness: cross-validation design and year",
          "Dashed: ADW. Random CV holds out gauges whose neighbours stay in training, which favours pure interpolation; "
          "spatial CV is the estimate for ungauged areas.", top=0.97)
    ladder_legend(fig, raw=False, loc="lower center", bbox_to_anchor=(0.36, 0.0))
    save(fig, "fusion_fig6_robustness")


# --------------------------------------------------------------------------------------------------
# Report tables: rounded strings, so the Typst file only lays them out.


def fmt(x, digits=3, signed=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:+.{digits}f}" if signed else f"{x:.{digits}f}"


def write_table(name: str, header: list[str], rows: list[list[str]]) -> None:
    """One report table as JSON records of pre-formatted strings.

    JSON rather than CSV because the repository ignores *.csv, and the report folder has to compile from a
    fresh checkout the way the other interim reports do.
    """
    path = REPORT_DIR / "tables" / f"{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    records = [dict(zip(header, row)) for row in rows]
    path.write_text(json.dumps(records, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"  wrote tables/{name}.json")


def report_tables(summary, stratified, ranking, fusion, anchor_metrics, comparison, detection, coefficients,
                  weights) -> None:
    methods = ["raw", "idw", "adw", "tps", "gwr", "auto", "residual_gwr", "mgwr", "blend_mgwr", "blend_agrenv_mgwr"]
    rows = []
    for anchor in ANCHORS:
        for method in methods:
            row = one(summary, scheme=SCHEME, product=anchor, method=method)
            if row is None:
                continue
            rows.append([anchor, METHOD_LABELS[method], fmt(row["RMSE"]), fmt(row["MAE"]), fmt(row["Bias"], signed=True),
                         fmt(row["r"]), fmt(row["CSI_0.1"]), fmt(row["CSI_2.5"]), fmt(row["CSI_8.0"]),
                         fmt(pct(row["RMSE_improvement_vs_raw"]), 1, True), fmt(pct(row["RMSE_improvement_vs_adw"]), 1, True)])
    write_table("overall_metrics", ["anchor", "method", "RMSE", "MAE", "Bias", "r", "CSI0.1", "CSI2.5", "CSI8",
                                    "vs_raw", "vs_adw"], rows)

    rows = []
    for stratum in ("overall", "moderate", "heavy"):
        for other in ("idw", "tps", "gwr"):
            if other == "gwr":
                row = one(ranking, method="gwr", baseline="adw", stratum=stratum)
                sign, base = 1.0, (row["RMSE_baseline"] if row else None)
                adw = row["RMSE_baseline"] if row else None
                rmse = row["RMSE_gwr"] if row else None
            else:
                row = one(ranking, method="adw", baseline=other, stratum=stratum)
                sign, base = -1.0, (row["RMSE_gwr"] if row else None)
                adw = row["RMSE_gwr"] if row else None
                rmse = row["RMSE_baseline"] if row else None
            if row is None:
                continue
            lo, hi = sorted([pct(sign * row["ci_low"] / base), pct(sign * row["ci_high"] / base)])
            rows.append([stratum, METHOD_LABELS[other], fmt(rmse), fmt(adw), fmt(pct(sign * row["delta_RMSE"] / base), 1, True),
                         f"({lo:+.1f}, {hi:+.1f})", fmt(row["pvalue_holm"], 3)])
    write_table("interpolation_ranking", ["stratum", "method", "RMSE", "RMSE_adw", "change", "ci", "p_holm"], rows)

    rows = []
    for stratum in ("overall", "heavy"):
        cells = [row for row in fusion if row["stratum"] == stratum]
        cells.sort(key=lambda r: -r["relative_improvement"])
        for row in cells:
            rows.append([stratum, row["product"], METHOD_LABELS[row["method"]], fmt(row["RMSE_method"]),
                         fmt(row["RMSE_adw"]), fmt(pct(row["relative_improvement"]), 1, True),
                         f"({pct(row['relative_ci_low']):+.1f}, {pct(row['relative_ci_high']):+.1f})",
                         fmt(row["pvalue_holm_family"], 3)])
    write_table("fusion_vs_adw", ["stratum", "anchor", "method", "RMSE", "RMSE_adw", "change", "ci", "p_holm_family"],
                rows)

    rows = []
    for anchor in ANCHORS + ["gauge_only"]:
        for method in ("raw", "auto", "mgwr", "blend_agrenv_mgwr", "adw"):
            overall = one(anchor_metrics, product=anchor, method=method, group="overall", level="all")
            if overall is None:
                continue
            csi8 = one(anchor_metrics, product=anchor, method=method, group="event_threshold", level="8.0")
            label = "Gauges only" if anchor == "gauge_only" else anchor
            rows.append([label, METHOD_LABELS[method], fmt(overall["RMSE"]), fmt(overall["Bias"], signed=True),
                         fmt(overall["r"]), fmt(csi8["CSI"] if csi8 else None)])
    write_table("anchor_metrics", ["anchor", "method", "RMSE", "Bias", "r", "CSI8"], rows)

    rows = [[row["product"], METHOD_LABELS[row["method"]], fmt(row["RMSE_product"]), fmt(row["RMSE_baseline"]),
             fmt(pct(row["relative_improvement"]), 1, True),
             f"({pct(row['relative_ci_low']):+.1f}, {pct(row['relative_ci_high']):+.1f})", fmt(row["pvalue_holm"], 3)]
            for row in comparison]
    write_table("anchor_vs_gpm", ["anchor", "method", "RMSE", "RMSE_gpm", "change", "ci", "p_holm"], rows)

    rows = []
    for anchor in ANCHORS:
        for method in ("raw", "tps", "auto", "blend_agrenv_mgwr"):
            for threshold in (2.5, 8.0, 16.0):
                cells = {m: one(detection, product=anchor, method=method, threshold=threshold, metric=m)
                         for m in ("POD", "FAR", "CSI")}
                if cells["CSI"] is None:
                    continue
                rows.append([anchor, METHOD_LABELS[method], f"{threshold:g}",
                             fmt(cells["POD"]["value_method"]), fmt(cells["FAR"]["value_method"]),
                             fmt(cells["CSI"]["value_method"]), fmt(cells["CSI"]["value_adw"]),
                             fmt(cells["CSI"]["delta"], 3, True),
                             f"({cells['CSI']['ci_low']:+.3f}, {cells['CSI']['ci_high']:+.3f})"])
    write_table("detection_skill", ["anchor", "method", "threshold", "POD", "FAR", "CSI", "CSI_adw", "delta", "ci"], rows)

    rows = [[row["scheme"], row["product"], fmt(row["intercept_mean"]), fmt(row["beta_fy4b_mean"]),
             fmt(row["beta_gpm_mean"]), fmt(row["beta_gsmap_mean"]),
             f"{row['beta_gpm_min']:.2f}–{row['beta_gpm_max']:.2f}", str(int(row["folds_fell_back"]))]
            for row in coefficients if row["product"] == "MERGED_OLS"]
    write_table("merge_weights", ["scheme", "product", "intercept", "FY4B", "GPM", "GSMaP", "GPM_range", "fell_back"],
                rows)

    # The agreement axis carries nine weights per fold; they stay in blend_weights.csv. The report shows the
    # single constant weight, which is the one a reader can interpret.
    rows = [[row["product"], fmt(row["lambda_mean"], 2), f"{row['lambda_min']:.1f}–{row['lambda_max']:.1f}",
             fmt(pct(row["mean_inner_improvement"]), 1, True)]
            for row in weights if row["method"] == "blend_mgwr"]
    write_table("blend_weights", ["anchor", "lambda_mean", "lambda_range", "inner_gain"], rows)

    rows = []
    for anchor in ANCHORS:
        for method in ("auto", "mgwr", "blend_mgwr", "blend_agrenv_mgwr"):
            for group, levels in (("rain_intensity", INTENSITY), ("nearest_train_km", DISTANCE)):
                for level, _ in levels:
                    change = value(stratified, "RMSE_improvement_vs_adw", product=anchor, method=method, group=group,
                                   level=level)
                    rows.append([anchor, METHOD_LABELS[method], group, level, fmt(pct(change), 1, True)])
    write_table("where_it_helps", ["anchor", "method", "group", "level", "change"], rows)

    rows = []
    for anchor in ANCHORS:
        for scheme in ("balanced_spatial", "random"):
            adw = value(summary, "RMSE", scheme=scheme, product=anchor, method="adw")
            cells = [scheme.replace("_", " "), anchor, fmt(adw)]
            for method in ("auto", "blend_agrenv_mgwr"):
                rmse = value(summary, "RMSE", scheme=scheme, product=anchor, method=method)
                cells += [fmt(rmse), fmt(pct((adw - rmse) / adw), 1, True)]
            years = [fmt(pct(value(stratified, "RMSE_improvement_vs_adw", product=anchor, method="blend_agrenv_mgwr",
                                   group="year", level=year)), 1, True) if scheme == "balanced_spatial" else "–"
                     for year in ("2022", "2023", "2024")]
            rows.append(cells + years)
    write_table("robustness", ["scheme", "anchor", "RMSE_adw", "RMSE_infold", "infold_vs_adw", "RMSE_blend",
                               "blend_vs_adw", "blend_2022", "blend_2023", "blend_2024"], rows)

    provenance = read_rows(RUN_DIR / "source_provenance.csv")
    source = (RUN_DIR / "source_run.txt").read_text(encoding="utf-8").strip()
    write_table("provenance", ["key", "value"],
                [["run_directory", Path(source).name]] + [[row["key"], str(row["value"])] for row in provenance])


def main() -> None:
    summary = read_rows(RUN_DIR / "method_summary.csv")
    stratified = read_rows(RUN_DIR / "stratified_summary.csv")
    fusion = read_rows(RUN_DIR / "fusion_vs_gauge_only.csv")
    ranking = read_rows(RUN_DIR / "interpolation_ranking.csv")
    anchor_metrics = read_rows(RUN_DIR / "anchor_metrics.csv")
    comparison = read_rows(RUN_DIR / "anchor_comparison.csv")
    detection = read_rows(RUN_DIR / "detection_skill_bootstrap.csv")
    coefficients = read_rows(RUN_DIR / "merged_anchor_coefficients.csv")
    weights = read_rows(RUN_DIR / "blend_weights.csv")

    fig1_benefit(summary)
    fig2_interpolators(ranking, summary)
    fig3_heatmap(fusion)
    fig4_where(stratified, detection)
    fig5_products(anchor_metrics, comparison, coefficients)
    fig6_robustness(summary, stratified)

    (REPORT_DIR / "figures").mkdir(parents=True, exist_ok=True)
    for path in sorted(FIG_DIR.glob("fusion_fig*.png")):
        shutil.copy2(path, REPORT_DIR / "figures" / path.name)
    report_tables(summary, stratified, ranking, fusion, anchor_metrics, comparison, detection, coefficients, weights)
    print(f"Copied figures and tables into {REPORT_DIR}")


if __name__ == "__main__":
    main()
