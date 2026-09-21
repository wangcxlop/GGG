"""Draw the no-rain evaluation figures and the report tables.

Reads the CSVs written by `scripts/run_no_rain_evaluation.jl` (output/no_rain_evaluation/) and
`scripts/run_no_rain_fusion_evaluation.jl` (its fusion/ subfolder), writes PNG (300 dpi) and PDF figures to
output/no_rain_evaluation/figures/, then copies the figures into Interim_results/figures/ and writes rounded
report tables to Interim_results/tables/norain_*.json, so the Typst report compiles from that folder alone:

    py -3.13 scripts/plot_no_rain_evaluation.py

Colour carries one job per figure. The products keep the temporal report's hues (FY4B blue, GPM orange,
GSMaP green); uncalibrated IMERG is GPM's orange, dashed or hollow - the same product before its gauge
adjustment - rather than a fourth hue, which would sit too close to FY4B's blue. The gauge is ink. The fusion
figures put methods on an axis and colour the five anchors with the validated slots 1-5 (the three products,
then amber and pink for the two merges). Magnitude maps use one blue ramp.
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
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "output" / "no_rain_evaluation"
FUSION_DIR = RUN_DIR / "fusion"
FIG_DIR = RUN_DIR / "figures"
REPORT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "Interim_results"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
PRODUCT_COLORS = {"FY4B": "#2a78d6", "GPM": "#eb6834", "GSMaP": "#1baf7a", "GPM_uncal": "#eb6834", "Gauge": INK}
PRODUCT_LABELS = {"FY4B": "FY4B", "GPM": "GPM", "GSMaP": "GSMaP", "GPM_uncal": "GPM uncalibrated", "Gauge": "Gauge"}
ANCHOR_COLORS = {"FY4B": "#2a78d6", "GPM": "#eb6834", "GSMaP": "#1baf7a", "MERGED_MEAN": "#eda100",
                 "MERGED_OLS": "#e87ba4"}
ANCHOR_LABELS = {"FY4B": "FY4B", "GPM": "GPM", "GSMaP": "GSMaP", "MERGED_MEAN": "Merged (equal)",
                 "MERGED_OLS": "Merged (OLS)"}
ANCHORS = list(ANCHOR_COLORS)
SEQUENTIAL_BLUES = ["#dcebfb", "#b7d4f6", "#86b6ef", "#5598e7", "#3987e5", "#256abf", "#1c5cab", "#104281"]
ORDINAL_BLUES = ["#86b6ef", "#5598e7", "#256abf", "#104281"]
SAMPLES = {"all_products": ["FY4B", "GPM", "GSMaP"], "gpm_gsmap_full": ["GPM", "GSMaP"],
           "gpm_cal_uncal": ["GPM", "GPM_uncal"]}
SAMPLE_TITLES = {"all_products": "Hours FY4B covers (three products)", "gpm_gsmap_full": "Full 2022–2024 record",
                 "gpm_cal_uncal": "Full record, calibrated vs uncalibrated IMERG"}
SEASONS = ["MAM", "JJA", "SON", "DJF"]
SCHEME = "balanced_spatial"
METHODS = ["raw", "adw", "idw", "tps", "gwr", "mgwr", "auto", "blend_mgwr", "blend_agrenv_mgwr"]
METHOD_LABELS = {
    "raw": "Raw satellite", "adw": "ADW (gauges only)", "idw": "IDW", "tps": "TPS", "gwr": "GWR (gauges only)",
    "auto": "GWR family, chosen in-fold", "residual_gwr": "Residual GWR", "mixed_gwr": "Mixed GWR", "mgwr": "MGWR",
    "blend_residual_gwr": "Residual GWR + ADW blend", "blend_mixed_gwr": "Mixed GWR + ADW blend",
    "blend_mgwr": "MGWR + ADW blend", "blend_agrenv_residual_gwr": "Residual GWR + agreement blend",
    "blend_agrenv_mixed_gwr": "Mixed GWR + agreement blend", "blend_agrenv_mgwr": "MGWR + agreement blend",
    "zero": "Always zero", "train_clim": "Training climatology", "hour_field_mean": "Hourly field mean",
}
STRING_COLUMNS = {"sample", "product", "stratifier", "level", "source", "pair", "scheme", "anchor", "method", "key",
                  "value", "station_id", "relief_region"}

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


# --------------------------------------------------------------------------------------------------
# Reading and small helpers


def read_rows(path: Path) -> list[dict]:
    """CSV rows with numeric-looking fields converted to float and true/false to bool."""
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            parsed = {}
            for key, value in row.items():
                if key in STRING_COLUMNS:
                    parsed[key] = value
                elif value in ("true", "false"):
                    parsed[key] = value == "true"
                else:
                    try:
                        parsed[key] = float(value)
                    except ValueError:
                        parsed[key] = value
            rows.append(parsed)
    return rows


def pick(rows: list[dict], **conditions) -> list[dict]:
    return [row for row in rows if all(row.get(k) == v for k, v in conditions.items())]


def one(rows: list[dict], **conditions) -> dict:
    hits = pick(rows, **conditions)
    if len(hits) != 1:
        raise KeyError(f"{len(hits)} rows match {conditions}")
    return hits[0]


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


def title(fig: plt.Figure, heading: str, note: str, left: float = 0.07) -> None:
    fig.suptitle(heading, x=left, ha="left", fontsize=9, y=0.975)
    fig.text(left, 0.925, note, ha="left", va="top", fontsize=6.3, color=INK_SECONDARY, linespacing=1.4)


def product_style(product: str) -> dict:
    """Line style of a product: uncalibrated IMERG is GPM's hue, dashed."""
    return {"color": PRODUCT_COLORS[product], "linestyle": "--" if product == "GPM_uncal" else "-", "linewidth": 1.5}


def product_legend(fig_or_ax, products, **kwargs) -> None:
    handles = [Line2D([], [], marker="o", markersize=4.5, markerfacecolor=SURFACE if p == "GPM_uncal" else PRODUCT_COLORS[p],
                      markeredgecolor=PRODUCT_COLORS[p], **{k: v for k, v in product_style(p).items() if k != "color"},
                      color=PRODUCT_COLORS[p], label=PRODUCT_LABELS[p]) for p in products]
    fig_or_ax.legend(handles=handles, frameon=False, fontsize=6.3, handlelength=2.2, **kwargs)


def dots(ax, y, values, color, lo=None, hi=None, hollow=False, size=24, horizontal=True):
    """Values as dots with optional interval whiskers; horizontal=True puts the value on x."""
    y, values = np.asarray(y, dtype=float), np.asarray(values, dtype=float)
    face = SURFACE if hollow else color
    if lo is not None:
        lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
        if horizontal:
            ax.hlines(y, lo, hi, color=color, linewidth=1.0, alpha=0.75, zorder=2)
        else:
            ax.vlines(y, lo, hi, color=color, linewidth=1.0, alpha=0.75, zorder=2)
    if horizontal:
        ax.scatter(values, y, s=size, facecolors=face, edgecolors=SURFACE if not hollow else color, linewidths=0.9, zorder=3)
    else:
        ax.scatter(y, values, s=size, facecolors=face, edgecolors=SURFACE if not hollow else color, linewidths=0.9, zorder=3)


def occ(rows, sample, product, stratifier, level, threshold) -> dict:
    return one(rows, sample=sample, product=product, stratifier=stratifier, level=level, threshold=threshold)


def pct(x: float, digits: int = 1) -> str:
    return "–" if not np.isfinite(x) else f"{100 * x:.{digits}f}"


def num(x: float, digits: int = 2) -> str:
    return "–" if not np.isfinite(x) else f"{x:.{digits}f}"


def interval(lo: float, hi: float, scale: float = 100, digits: int = 1) -> str:
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return ""
    return f"({scale * lo:.{digits}f}, {scale * hi:.{digits}f})"


def write_table(name: str, rows: list[dict]) -> None:
    """A report table as a JSON list of records; `*.csv` is gitignored, JSON is not."""
    path = REPORT_DIR / "tables" / f"norain_{name}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"  wrote tables/norain_{name}.json ({len(rows)} rows)")


# --------------------------------------------------------------------------------------------------
# Figure 1: the no-rain record itself


def fig1_anatomy(occurrence: list[dict]) -> None:
    sample, product = "gpm_gsmap_full", "GPM"
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 3.0), gridspec_kw={"width_ratios": [1.1, 1.5, 0.9, 0.9]})
    fig.subplots_adjust(left=0.06, right=0.99, top=0.74, bottom=0.2, wspace=0.62)

    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    share = [occ(occurrence, sample, product, "month", m, 0.1)["dry_share"] for m in months]
    ax = axes[0]
    clean_axis(ax)
    ax.bar(range(12), share, color=ORDINAL_BLUES[1], width=0.72, zorder=2)
    ax.set_xticks(range(12))
    ax.set_xticklabels([m[0] for m in months], fontsize=6)
    ax.set_ylim(0.8, 1.0)
    ax.set_title("Dry share by month", loc="left")

    ax = axes[1]
    clean_axis(ax)
    hours = [f"{h:02d}" for h in range(24)]
    for s_index, season in enumerate(SEASONS):
        values = [occ(occurrence, sample, product, "season_hour", f"{season}|{h}", 0.1)["dry_share"] for h in hours]
        ax.plot(range(24), values, color=ORDINAL_BLUES[s_index], linewidth=1.4, label=season, zorder=2)
    ax.set_xticks(range(0, 24, 3))
    ax.set_xlabel("Hour of day (BJT, hour-ending)")
    ax.set_title("Dry share by hour of day and season", loc="left")
    ax.legend(frameon=False, fontsize=6, ncol=4, loc="lower left")

    for ax, stratifier, heading in ((axes[2], "rain_proximity", "Dry hours by distance\nto the station's rain"),
                                    (axes[3], "network", "Dry hours by the\nnetwork's state")):
        clean_axis(ax, "x")
        rows = [r for r in pick(occurrence, sample=sample, product=product, stratifier=stratifier, threshold=0.1)]
        total = sum(r["n_dry"] for r in rows)
        levels = [r["level"] for r in rows]
        values = [r["n_dry"] / total for r in rows]
        y = np.arange(len(levels))[::-1]
        ax.barh(y, values, color=ORDINAL_BLUES[2], height=0.66, zorder=2)
        for yy, v in zip(y, values):
            ax.text(v + 0.01, yy, f"{100 * v:.0f}%", va="center", fontsize=6, color=INK_SECONDARY)
        ax.set_yticks(y)
        ax.set_yticklabels(levels, fontsize=6)
        ax.set_xlim(0, max(values) * 1.35)
        ax.set_xticks([])
        ax.set_title(heading, loc="left")
    title(fig, "The no-rain record: when and where the gauges are dry",
          "Full 2022–2024 record, 237 gauges; a dry hour is gauge < 0.1 mm/h (in practice a reading of 0.0)."
          "Distances to rain are to the same gauge's nearest wet hour;\nhours next to a missing gauge value are left out. "
          "Network state: share of reporting gauges that are wet in that hour (hours with < 90% reporting left out).",
          left=0.06)
    save(fig, "norain_fig1_anatomy")


# --------------------------------------------------------------------------------------------------
# Figure 2: exceedance on gauge-dry hours


def fig2_exceedance(exceedance: list[dict], baseline: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.76, bottom=0.16, wspace=0.08)
    for ax, sample, products in ((axes[0], "all_products", ["FY4B", "GPM", "GSMaP"]),
                                 (axes[1], "gpm_gsmap_full", ["GPM", "GSMaP", "GPM_uncal"])):
        clean_axis(ax, "both")
        for product in products:
            source = "gpm_cal_uncal" if product == "GPM_uncal" else sample
            rows = sorted(pick(exceedance, sample=source, product=product, stratifier="all", level="all"),
                          key=lambda r: r["threshold"])
            ax.plot([r["threshold"] for r in rows], [r["POFD"] for r in rows], marker="o", markersize=3,
                    markerfacecolor=SURFACE if product == "GPM_uncal" else PRODUCT_COLORS[product],
                    markeredgecolor=PRODUCT_COLORS[product], **product_style(product), zorder=3)
        near = one(baseline, sample=sample, level="0-5 km")
        ax.axhline(near["POFD"], color=MUTED, linewidth=0.9, linestyle=":", zorder=1)
        ax.text(0.011, near["POFD"] * 1.12, f"a gauge 0–5 km away: {100 * near['POFD']:.1f}%", fontsize=6, color=INK_SECONDARY)
        ax.axvline(0.5, color=AXIS, linewidth=0.8, zorder=1)
        ax.text(0.52, 0.2, "stricter\nthreshold", fontsize=5.8, color=MUTED, va="top")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Estimate threshold t (mm/h)")
        ax.set_title(SAMPLE_TITLES[sample], loc="left")
        ax.set_xticks([0.01, 0.1, 0.5, 1, 4])
        ax.set_xticklabels(["0.01", "0.1", "0.5", "1", "4"])
        product_legend(ax, products, loc="lower left")
    axes[0].set_ylabel("P(estimate ≥ t | gauge dry)")
    axes[0].yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    title(fig, "How often each product reports rain when the gauge is dry",
          "Share of gauge-dry station-hours on which the product reaches each threshold (all seasons, all gauges). "
          "Dotted: the same rate for a neighbouring gauge\n0–5 km away used as the estimate - the floor a point-to-area "
          "comparison cannot get under.", left=0.08)
    save(fig, "norain_fig2_exceedance")


# --------------------------------------------------------------------------------------------------
# Figure 3: the contingency scores

CONTINGENCY = [("POFD", "False-alarm rate\nP(wet | gauge dry)", 0.0, True),
               ("FAR", "False-alarm ratio\nP(gauge dry | wet)", 0.0, True),
               ("NPV", "Dry reliability\nP(gauge dry | dry)", 1.0, True),
               ("freq_bias", "Frequency bias\nwet hours / gauge wet hours", 1.0, True),
               ("HSS", "Heidke skill score", 1.0, True)]


def fig3_contingency(occurrence: list[dict]) -> None:
    fig, axes = plt.subplots(1, len(CONTINGENCY), figsize=(10.0, 2.9), sharey=True)
    fig.subplots_adjust(left=0.1, right=0.99, top=0.7, bottom=0.14, wspace=0.18)
    entries = [("all_products", p) for p in SAMPLES["all_products"]] + \
              [("gpm_gsmap_full", p) for p in SAMPLES["gpm_gsmap_full"]] + [("gpm_cal_uncal", "GPM_uncal")]
    labels = [f"{PRODUCT_LABELS[p]}" + (" · FY4B hours" if s == "all_products" else "") for s, p in entries]
    y = np.arange(len(entries))[::-1].astype(float)
    y[:3] += 0.4
    for ax, (column, heading, ideal, with_ci) in zip(axes, CONTINGENCY):
        clean_axis(ax, "x")
        ax.axvline(ideal, color=AXIS, linewidth=0.9, zorder=1)
        for k, (sample, product) in enumerate(entries):
            for threshold, hollow in ((0.1, False), (0.5, True)):
                r = occ(occurrence, sample, product, "all", "all", threshold)
                dots(ax, [y[k] + (0.13 if threshold == 0.5 else -0.13)], [r[column]], PRODUCT_COLORS[product],
                     [r[f"{column}_lo"]], [r[f"{column}_hi"]], hollow=hollow, size=20)
        ax.set_title(heading, loc="left", fontsize=6.8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=6.3)
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=INK_SECONDARY,
                      markeredgecolor=INK_SECONDARY, label="estimate ≥ 0.1 mm/h"),
               Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=SURFACE,
                      markeredgecolor=INK_SECONDARY, label="estimate ≥ 0.5 mm/h")]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.99), frameon=False, fontsize=6.3, ncol=2)
    title(fig, "Occurrence scores: the dry side of the contingency table",
          "Gauge wet at ≥ 0.1 mm/h. Whiskers: 95% day-block bootstrap. "
          "Grey line: a perfect score.", left=0.1)
    save(fig, "norain_fig3_contingency")


# --------------------------------------------------------------------------------------------------
# Figure 4: diurnal cycle of false alarms


def fig4_diurnal(occurrence: list[dict]) -> None:
    fig, axes = plt.subplots(1, 5, figsize=(10.0, 2.9), sharey=True)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.72, bottom=0.17, wspace=0.08)
    hours = [f"{h:02d}" for h in range(24)]
    ax = axes[0]
    clean_axis(ax)
    for product in SAMPLES["all_products"]:
        values = [occ(occurrence, "all_products", product, "hour", h, 0.5)["POFD"] for h in hours]
        ax.plot(range(24), values, **product_style(product), zorder=2)
    ax.set_title("All seasons · FY4B hours", loc="left")
    ax.set_ylabel("False-alarm rate, estimate ≥ 0.5 mm/h")
    for ax, season in zip(axes[1:], SEASONS):
        clean_axis(ax)
        for product in SAMPLES["gpm_gsmap_full"]:
            values = [occ(occurrence, "gpm_gsmap_full", product, "season_hour", f"{season}|{h}", 0.5)["POFD"] for h in hours]
            ax.plot(range(24), values, **product_style(product), zorder=2)
        ax.set_title(f"{season} · full record", loc="left")
    for ax in axes:
        ax.set_xticks(range(0, 24, 6))
        ax.set_xlim(-0.5, 23.5)
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    axes[2].set_xlabel("Hour of day (BJT, hour-ending)")
    product_legend(fig, ["FY4B", "GPM", "GSMaP"], loc="upper right", bbox_to_anchor=(0.99, 0.99), ncol=3)
    title(fig, "When in the day the false alarms happen",
          "Share of gauge-dry station-hours on which the product reports ≥ 0.5 mm/h, by hour of day. "
          "FY4B never provides the 23:00 label.", left=0.06)
    save(fig, "norain_fig4_diurnal")


# --------------------------------------------------------------------------------------------------
# Figure 5: station maps


def fig5_maps(occurrence: list[dict], stations: list[dict]) -> None:
    sample = "all_products"
    products = SAMPLES[sample]
    fig, axes = plt.subplots(2, 3, figsize=(8.4, 6.0))
    fig.subplots_adjust(left=0.06, right=0.9, top=0.86, bottom=0.05, wspace=0.08, hspace=0.18)
    lon = {s["station_id"]: s["lon"] for s in stations}
    lat = {s["station_id"]: s["lat"] for s in stations}
    # Eight bounds: seven bins plus the "max" extension use the ramp's eight colours.
    # The volume share, not mm per year: the FY4B-hours sample is weighted to summer.
    rows_spec = [("POFD", 0.5, "False-alarm rate ≥ 0.5 mm/h (%)", [0, 2.5, 3, 3.5, 4, 4.5, 5, 6], 100),
                 ("dry_volume_share", 0.1, "Share of the product's rain on gauge-dry hours (%)",
                  [0, 25, 30, 35, 40, 50, 60, 70], 100)]
    cmap = ListedColormap(SEQUENTIAL_BLUES)
    for i, (column, threshold, label, bounds, scale) in enumerate(rows_spec):
        norm = BoundaryNorm(bounds, cmap.N, extend="max")
        for j, product in enumerate(products):
            ax = axes[i, j]
            rows = pick(occurrence, sample=sample, product=product, stratifier="station", threshold=threshold)
            x = [lon[r["level"]] for r in rows]
            yy = [lat[r["level"]] for r in rows]
            values = [scale * r[column] for r in rows]
            ax.scatter(x, yy, c=values, cmap=cmap, norm=norm, s=16, edgecolors=INK_SECONDARY, linewidths=0.3, zorder=2)
            ax.set_aspect(1 / np.cos(np.radians(32.2)), adjustable="box")
            ax.tick_params(length=2, pad=1.5, labelsize=5.6)
            ax.set_xticks([110, 110.5, 111])
            ax.set_yticks([31.5, 32, 32.5, 33])
            if i == 0:
                ax.set_title(PRODUCT_LABELS[product], loc="left")
            if j > 0:
                ax.set_yticklabels([])
            med = np.median(values)
            ax.text(0.02, 0.03, f"median {med:.1f}%", transform=ax.transAxes, fontsize=6, color=INK_SECONDARY)
        cax = fig.add_axes([0.915, 0.53 - 0.44 * i, 0.012, 0.3])
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, extend="max")
        bar.outline.set_edgecolor(AXIS)
        bar.ax.tick_params(length=2, labelsize=5.8)
        bar.set_label(label, fontsize=6.3, color=INK_SECONDARY)
    title(fig, "Where each product rains on gauge-dry hours",
          "Hours FY4B covers, all seasons; each dot is a gauge and the pixel over it. Bottom: the share of the "
          "product's total at that pixel\nthat falls on hours the gauge is dry.", left=0.06)
    save(fig, "norain_fig5_station_maps")


# --------------------------------------------------------------------------------------------------
# Figure 6: is the false rain really false?

CONTEXT_PANELS = [("rain_proximity", "Distance to the station's own rain"),
                  ("rain_side", "Side of the station's rain (± 3 h)"),
                  ("neighbours", "Gauges within 15 km"),
                  ("network", "Whole network")]


def fig6_context(occurrence: list[dict]) -> None:
    sample = "all_products"
    products = SAMPLES[sample]
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 3.1), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.72, bottom=0.15, wspace=0.95)
    for ax, (stratifier, heading) in zip(axes, CONTEXT_PANELS):
        clean_axis(ax, "x")
        levels = [r["level"] for r in pick(occurrence, sample=sample, product="GPM", stratifier=stratifier, threshold=0.5)]
        y = np.arange(len(levels))[::-1].astype(float)
        for k, product in enumerate(products):
            offset = (k - 1) * 0.22
            rows = [occ(occurrence, sample, product, stratifier, level, 0.5) for level in levels]
            dots(ax, y + offset, [r["POFD"] for r in rows], PRODUCT_COLORS[product],
                 [r["POFD_lo"] for r in rows], [r["POFD_hi"] for r in rows], size=18)
        shares = [occ(occurrence, sample, "GPM", stratifier, level, 0.5)["n_dry"] for level in levels]
        total = sum(shares)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{lv}  ({100 * s / total:.0f}%)" for lv, s in zip(levels, shares)], fontsize=6)
        ax.set_title(heading, loc="left")
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    axes[1].set_xlabel("False-alarm rate, estimate ≥ 0.5 mm/h (log scale)", x=1.2)
    product_legend(fig, products, loc="upper right", bbox_to_anchor=(0.99, 0.99), ncol=3)
    title(fig, "Is the satellite's rain on dry gauge hours really false?",
          "Hours FY4B covers. Gauge-dry hours grouped by what the gauges show around them; in brackets, each group's "
          "share of the dry hours. Whiskers: 95% day-block bootstrap.", left=0.08)
    save(fig, "norain_fig6_context")


# --------------------------------------------------------------------------------------------------
# Figure 7: the gauge-to-gauge floor and meteorology


def fig7_floor_and_air(occurrence: list[dict], baseline: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.1), gridspec_kw={"width_ratios": [1.2, 1.0, 1.0]})
    fig.subplots_adjust(left=0.07, right=0.99, top=0.72, bottom=0.2, wspace=0.3)
    ax = axes[0]
    clean_axis(ax)
    sample = "all_products"
    bins = [r for r in pick(baseline, sample=sample) if r["level"] != "nearest gauge"]
    x = np.arange(len(bins))
    dots(ax, x, [r["POFD"] for r in bins], INK, [r["POFD_lo"] for r in bins], [r["POFD_hi"] for r in bins],
         size=22, horizontal=False)
    # Reference lines stop short of their labels; labels that would overprint are pushed apart.
    levels = sorted((occ(occurrence, sample, p, "all", "all", 0.5)["POFD"], p) for p in SAMPLES[sample])
    label_y = [v for v, _ in levels]
    for k in range(1, len(label_y)):
        label_y[k] = max(label_y[k], label_y[k - 1] + 0.0017)
    for (value, product), y_text in zip(levels, label_y):
        ax.hlines(value, -0.5, len(bins) - 0.55, color=PRODUCT_COLORS[product], linewidth=1.2, zorder=1)
        ax.text(len(bins) - 0.5, y_text, f"{product} {100 * value:.1f}%", color=PRODUCT_COLORS[product],
                fontsize=6, va="center")
    ax.set_xticks(x)
    ax.set_xticklabels([r["level"] for r in bins], fontsize=6)
    ax.set_xlim(-0.5, len(bins) + 0.7)
    ax.set_ylim(0, None)
    ax.set_xlabel("Separation of the two gauges")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    ax.set_title("A neighbouring gauge as the estimate", loc="left")
    ax.set_ylabel("False-alarm rate")

    for ax, stratifier, heading, xlabel in ((axes[1], "dewpoint_depression", "Air near saturation?", "2 m dew-point depression (°C)"),
                                            (axes[2], "t2m", "Cold or warm?", "2 m temperature (°C)")):
        clean_axis(ax)
        levels = [r["level"] for r in pick(occurrence, sample=sample, product="GPM", stratifier=stratifier, threshold=0.5)]
        x = np.arange(len(levels))
        for k, product in enumerate(SAMPLES[sample]):
            rows = [occ(occurrence, sample, product, stratifier, level, 0.5) for level in levels]
            ax.plot(x, [r["POFD"] for r in rows], marker="o", markersize=3.5, **product_style(product), zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels(levels, fontsize=6)
        ax.set_xlabel(xlabel)
        ax.set_title(heading, loc="left")
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
        ax.set_ylim(0, None)
    product_legend(axes[2], SAMPLES[sample], loc="upper left")
    title(fig, "Two physical baselines: gauges next to each other, and the air at the gauge",
          "Hours FY4B covers. Left: how often a gauge reports rain when a neighbour at the given separation is dry, "
          "against each product's rate at ≥ 0.5 mm/h (lines).\nMiddle and right: the products' false-alarm rate by "
          "ERA5-Land 2 m dew-point depression and temperature at the gauge.", left=0.07)
    save(fig, "norain_fig7_floor_meteorology")


# --------------------------------------------------------------------------------------------------
# Figure 8: agreement between products, and the network-dry hours


def fig8_agreement(occurrence: list[dict], extent: list[dict]) -> None:
    sample = "all_products"
    products = SAMPLES[sample]
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.0))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.72, bottom=0.18, wspace=0.35)
    levels = ["0 others wet", "1 other wet", "2 others wet"]
    # The joint distribution of "how many products are wet", recovered from FY4B's conditional rows.
    r = [occ(occurrence, sample, "FY4B", "other_products", lv, 0.1) for lv in levels]
    n = [x["n_dry"] for x in r]
    p = [x["POFD"] for x in r]
    joint = [n[0] * (1 - p[0]), n[0] * p[0] + n[1] * (1 - p[1]), n[1] * p[1] + n[2] * (1 - p[2]), n[2] * p[2]]
    total = sum(joint)
    ax = axes[0]
    clean_axis(ax)
    ax.bar(range(4), [j / total for j in joint], color=ORDINAL_BLUES[2], width=0.66, zorder=2)
    for k, j in enumerate(joint):
        ax.text(k, j / total, f"{100 * j / total:.1f}%", ha="center", va="bottom", fontsize=6, color=INK_SECONDARY)
    ax.set_yscale("log")
    ax.set_xticks(range(4))
    ax.set_xticklabels(["none", "one", "two", "all three"])
    ax.set_xlabel("Products reporting ≥ 0.1 mm/h")
    ax.set_title("Gauge-dry hours by product agreement", loc="left")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))

    ax = axes[1]
    clean_axis(ax)
    for product in products:
        rows = [occ(occurrence, sample, product, "other_products", lv, 0.1) for lv in levels]
        ax.plot(range(3), [x["POFD"] for x in rows], marker="o", markersize=3.5, **product_style(product), zorder=2)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["none", "one", "both"])
    ax.set_xlabel("Of the other two products, wet")
    ax.set_title("P(product wet | gauge dry, others' state)", loc="left")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    product_legend(ax, products, loc="upper left")

    ax = axes[2]
    clean_axis(ax, "x")
    cut = [("share_hours_gt_0", "any gauge pixel"), ("share_hours_gt_5", "> 5% of pixels"),
           ("share_hours_gt_25", "> 25% of pixels"), ("share_hours_gt_50", "> 50% of pixels")]
    y = np.arange(len(cut))[::-1].astype(float)
    for k, product in enumerate(products):
        row = one(extent, sample=sample, product=product)
        ax.barh(y + (k - 1) * 0.26, [row[c] for c, _ in cut], height=0.24, color=PRODUCT_COLORS[product], zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label in cut], fontsize=6)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    ax.set_xlabel("Share of network-dry hours")
    n_hours = int(one(extent, sample=sample, product="GPM")["n_hours"])
    ax.set_title(f"Rain on hours no gauge is wet ({n_hours:,} h)", loc="left")
    title(fig, "Do the products agree when they rain on a dry gauge?",
          "Hours FY4B covers. Left, middle: gauge-dry station-hours at ≥ 0.1 mm/h. Right: hours on which all reporting "
          "gauges are dry (≥ 90% reporting), and how widely each product rains.", left=0.07)
    save(fig, "norain_fig8_agreement")


# --------------------------------------------------------------------------------------------------
# Figure 9: dry spells and dry days


def survival(hist: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    """P(spell length >= lower bound) from binned complete spells, weighted by the hours they hold."""
    rows = sorted(hist, key=lambda r: r["bin_lower"])
    cells = np.array([r["cells_in_spells"] for r in rows], dtype=float)
    lower = np.array([r["bin_lower"] for r in rows], dtype=float)
    tail = cells[::-1].cumsum()[::-1] / cells.sum()
    return lower, tail


def fig9_spells(hourly_hist: list[dict], daily: list[dict], daily_summary: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.2), gridspec_kw={"width_ratios": [1.2, 1.0, 1.0]})
    fig.subplots_adjust(left=0.07, right=0.99, top=0.72, bottom=0.17, wspace=0.32)
    ax = axes[0]
    clean_axis(ax, "both")
    sample = "gpm_gsmap_full"
    for source in ["Gauge", "GPM", "GSMaP"]:
        lower, tail = survival(pick(hourly_hist, sample=sample, source=source, threshold=0.1))
        style = product_style(source) if source != "Gauge" else {"color": INK, "linewidth": 1.8}
        ax.step(lower, tail, where="post", **style, label=PRODUCT_LABELS[source], zorder=3)
    lower, tail = survival(pick(hourly_hist, sample="gpm_cal_uncal", source="GPM_uncal", threshold=0.1))
    ax.step(lower, tail, where="post", **product_style("GPM_uncal"), label="GPM uncalibrated", zorder=3)
    ax.set_xscale("log")
    ax.set_xticks([1, 3, 6, 12, 24, 72, 168, 384])
    ax.set_xticklabels(["1", "3", "6", "12", "24", "72", "168", "384"])
    ax.set_xlabel("Dry spell length L (hours, estimate < 0.1 mm/h)")
    ax.set_ylabel("Share of dry hours in spells ≥ L")
    ax.set_title("Hourly dry spells · full record", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="lower left")

    ax = axes[1]
    clean_axis(ax, "x")
    entries = [("all_products", p) for p in SAMPLES["all_products"]] + [("gpm_gsmap_full", p) for p in SAMPLES["gpm_gsmap_full"]]
    y = np.arange(len(entries))[::-1].astype(float)
    for k, (s, product) in enumerate(entries):
        for threshold, hollow, offset in ((0.1, False, -0.14), (1.0, True, 0.14)):
            r = one(daily, sample=s, product=product, stratifier="all", level="all", threshold=threshold)
            dots(ax, [y[k] + offset], [r["freq_bias"]], PRODUCT_COLORS[product], [r["freq_bias_lo"]], [r["freq_bias_hi"]],
                 hollow=hollow, size=20)
    ax.axvline(1.0, color=AXIS, linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([PRODUCT_LABELS[p] + (" · FY4B hours" if s == "all_products" else "") for s, p in entries], fontsize=6)
    ax.set_title("Wet-day frequency bias", loc="left")
    ax.set_xlabel("Product wet days / gauge wet days")
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=INK_SECONDARY,
                      markeredgecolor=INK_SECONDARY, label="wet day ≥ 0.1 mm"),
               Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=SURFACE,
                      markeredgecolor=INK_SECONDARY, label="wet day ≥ 1 mm")]
    ax.legend(handles=handles, frameon=False, fontsize=6, loc="upper right")

    ax = axes[2]
    clean_axis(ax, "both")
    sample = "gpm_gsmap_full"
    for product in SAMPLES[sample]:
        gauge = {r["level"]: r["max_length"] for r in pick(daily_summary, sample=sample, pair=product, source="Gauge", threshold=1.0)
                 if r["level"] != "all"}
        est = {r["level"]: r["max_length"] for r in pick(daily_summary, sample=sample, pair=product, source=product, threshold=1.0)
               if r["level"] != "all"}
        keys = sorted(gauge)
        ax.scatter([gauge[k] for k in keys], [est[k] for k in keys], s=9, color=PRODUCT_COLORS[product], alpha=0.8,
                   edgecolors="none", zorder=2, label=PRODUCT_LABELS[product])
    limit = ax.get_xlim()[1]
    ax.plot([0, 200], [0, 200], color=AXIS, linewidth=0.9, zorder=1)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.set_xlabel("Gauge longest dry spell (days < 1 mm)")
    ax.set_ylabel("Product longest dry spell (days)")
    ax.set_title("Longest dry spell at each gauge · full record", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    title(fig, "Dry spells and dry days",
          "Left: complete hourly dry spells (spells cut by missing hours are not counted), full record. Middle: 08–08 "
          "BJT days; FY4B-hour days count when ≥ 20 of 24 hours are scored,\nsummed over the same hours for gauge and "
          "product. Right: longest run of days < 1 mm at each gauge over 2022–2024.", left=0.07)
    save(fig, "norain_fig9_spells")


# --------------------------------------------------------------------------------------------------
# Figure 10: calibrated vs uncalibrated IMERG


def fig10_calibration(occurrence: list[dict], exceedance: list[dict], lags: list[dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.0), gridspec_kw={"width_ratios": [1.2, 1.2, 0.8]})
    fig.subplots_adjust(left=0.07, right=0.99, top=0.72, bottom=0.18, wspace=0.32)
    sample = "gpm_cal_uncal"
    ax = axes[0]
    clean_axis(ax, "both")
    for product in SAMPLES[sample]:
        rows = sorted(pick(exceedance, sample=sample, product=product, stratifier="all", level="all"), key=lambda r: r["threshold"])
        ax.plot([r["threshold"] for r in rows], [r["POFD"] for r in rows], marker="o", markersize=3,
                markerfacecolor=SURFACE if product == "GPM_uncal" else PRODUCT_COLORS[product],
                markeredgecolor=PRODUCT_COLORS[product], **product_style(product), zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Threshold t (mm/h)")
    ax.set_xticks([0.01, 0.1, 0.5, 1, 4])
    ax.set_xticklabels(["0.01", "0.1", "0.5", "1", "4"])
    ax.set_ylabel("P(estimate ≥ t | gauge dry)")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    ax.set_title("Exceedance on gauge-dry hours", loc="left")
    product_legend(ax, SAMPLES[sample], loc="lower left")

    ax = axes[1]
    clean_axis(ax)
    x = np.arange(len(SEASONS) + 1)
    for k, product in enumerate(SAMPLES[sample]):
        for threshold, offset in ((0.1, -0.18), (0.5, 0.18)):
            rows = [occ(occurrence, sample, product, "all", "all", threshold)] + \
                   [occ(occurrence, sample, product, "season", s, threshold) for s in SEASONS]
            dots(ax, x + offset + (k - 0.5) * 0.14, [r["POFD"] for r in rows], PRODUCT_COLORS[product],
                 [r["POFD_lo"] for r in rows], [r["POFD_hi"] for r in rows],
                 hollow=product == "GPM_uncal", size=18, horizontal=False)
    ax.set_xticks(x)
    ax.set_xticklabels(["All"] + SEASONS)
    ax.set_ylim(0, None)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    ax.set_title("False-alarm rate by season (left pair ≥ 0.1, right pair ≥ 0.5)", loc="left")

    ax = axes[2]
    clean_axis(ax)
    ax.bar([r["lag"] for r in lags], [r["r"] for r in lags], color=ORDINAL_BLUES[1], width=0.7, zorder=2)
    ax.set_xlabel("Lag of uncalibrated series (h)")
    ax.set_ylabel("Pooled correlation r (unitless)")
    ax.set_title("Clock check: pooled r", loc="left")
    title(fig, "Does IMERG's gauge calibration remove the dry-hour rain?",
          "Full record, cells where both IMERG series and the gauge report. The calibrated series is the Final run's "
          "`precipitation`; the uncalibrated one is its\n`precipitationUncal`, downloaded separately and placed on the "
          "gauge clock (UTC hour start + 9 h).", left=0.07)
    save(fig, "norain_fig10_calibration")


# --------------------------------------------------------------------------------------------------
# Figure 13: silent gauges


def fig13_silent(silent: list[dict], occurrence: list[dict]) -> None:
    flagged = [r for r in silent if r["flagged"]]
    fig = plt.figure(figsize=(10.0, 3.6))
    ax_time = fig.add_axes([0.07, 0.17, 0.46, 0.58])
    ax_rate = fig.add_axes([0.66, 0.17, 0.32, 0.58])
    clean_axis(ax_time, "x")
    first = {}
    for r in flagged:
        first[r["station_id"]] = min(first.get(r["station_id"], 1e9), r["first_day"])
    order = sorted(first, key=lambda sid: first[sid])
    row_of = {sid: k for k, sid in enumerate(order)}
    for r in flagged:
        y = len(order) - 1 - row_of[r["station_id"]]
        ax_time.barh(y, r["last_day"] - r["first_day"] + 1, left=r["first_day"], height=0.7, color=ORDINAL_BLUES[2], zorder=2)
    year_starts = [1, 366, 731, 1097]
    ax_time.set_xticks([183, 548, 914])
    ax_time.set_xticklabels(["2022", "2023", "2024"])
    for x in year_starts:
        ax_time.axvline(x, color=AXIS, linewidth=0.6, zorder=1)
    ax_time.set_yticks([])
    ax_time.set_ylabel(f"{len(order)} gauges with a flagged spell")
    ax_time.set_xlim(1, 1097)
    ax_time.set_title("When the flagged gauges were silent (08-08 days)", loc="left")

    clean_axis(ax_rate, "x")
    entries = [("all_products", p) for p in SAMPLES["all_products"]] + [("gpm_gsmap_full", p) for p in SAMPLES["gpm_gsmap_full"]]
    y = np.arange(len(entries))[::-1].astype(float)
    for k, (sample, product) in enumerate(entries):
        for level, hollow, offset in (("gauge working", False, -0.14), ("silent-gauge spell", True, 0.14)):
            r = occ(occurrence, sample, product, "gauge_screen", level, 0.5)
            dots(ax_rate, [y[k] + offset], [r["POFD"]], PRODUCT_COLORS[product], [r["POFD_lo"]], [r["POFD_hi"]],
                 hollow=hollow, size=20)
    ax_rate.set_xscale("log")
    # The span is under one decade, so the labelled ticks are log-minor ones.
    ax_rate.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    ax_rate.xaxis.set_minor_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    ax_rate.set_yticks(y)
    ax_rate.set_yticklabels([PRODUCT_LABELS[p] + (" · FY4B hours" if s == "all_products" else "") for s, p in entries],
                            fontsize=6.3)
    ax_rate.set_xlabel("False-alarm rate, estimate ≥ 0.5 mm/h (log scale)")
    ax_rate.set_title("False alarms at working and silent gauges", loc="left")
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=INK_SECONDARY,
                      markeredgecolor=INK_SECONDARY, label="gauge working"),
               Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=SURFACE,
                      markeredgecolor=INK_SECONDARY, label="inside a silent-gauge spell")]
    ax_rate.legend(handles=handles, frameon=False, fontsize=6, loc="upper left")
    title(fig, "Some gauge zeros are not dry hours: silent gauges",
          "A run of ≥ 7 days < 1 mm is flagged when the median of the gauge's 3 nearest neighbours reaches 5 mm on "
          "≥ 3 of its days (a working gauge stays < 1 mm on only a few percent\nof such days). Right: the products' "
          "false-alarm rate on dry hours inside and outside the flagged spells; whiskers: 95% day-block bootstrap.",
          left=0.07)
    save(fig, "norain_fig13_silent_gauges")


# --------------------------------------------------------------------------------------------------
# Figure 11: the benchmark's methods on dry hours


def fig11_fusion(summary: list[dict], paired: list[dict]) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 3.6), sharey=True)
    fig.subplots_adjust(left=0.13, right=0.99, top=0.7, bottom=0.14, wspace=0.14)
    y = np.arange(len(METHODS))[::-1].astype(float)
    panels = [("spurious_mm_per_year", "Rain on dry hours (mm/yr)", None),
              ("POFD", "Dry hours predicted ≥ 0.1", 0.1),
              ("zero_share_dry", "Dry hours predicted exactly 0", None),
              ("dry_sse_share", "Dry hours' share of squared error", None)]
    for ax, (column, heading, _) in zip(axes, panels):
        clean_axis(ax, "x")
        for k, anchor in enumerate(ANCHORS):
            offset = (k - 2) * 0.13
            for m_index, method in enumerate(METHODS):
                if method in ("adw", "idw", "tps", "gwr") and anchor != "GPM":
                    continue
                r = one(summary, scheme=SCHEME, anchor=anchor, method=method, stratifier="all", threshold=0.1)
                dy = 0 if method in ("adw", "idw", "tps", "gwr") else offset
                color = INK if method in ("adw", "idw", "tps", "gwr") else ANCHOR_COLORS[anchor]
                lo = [r.get(f"{column}_lo", np.nan)] if f"{column}_lo" in r else None
                hi = [r.get(f"{column}_hi", np.nan)] if f"{column}_hi" in r else None
                dots(ax, [y[m_index] + dy], [r[column]], color, lo, hi, size=14)
        ax.set_title(heading, loc="left", fontsize=6.8)
        if column != "spurious_mm_per_year":
            ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{100 * v:g}%"))
    gauge_rate = one(summary, scheme=SCHEME, anchor="GPM", method="raw", stratifier="all", threshold=0.1)["gauge_mm_per_year"]
    axes[0].set_xlabel(f"mm per year (gauge total: {gauge_rate:.0f})", fontsize=6)
    axes[1].set_xlabel("% of dry hours", fontsize=6)
    axes[2].set_xlabel("% of dry hours", fontsize=6)
    axes[3].set_xlabel("% of total squared error", fontsize=6)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([METHOD_LABELS[m] for m in METHODS], fontsize=6.3)
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=ANCHOR_COLORS[a],
                      markeredgecolor=SURFACE, label=ANCHOR_LABELS[a]) for a in ANCHORS]
    handles.append(Line2D([], [], marker="o", linestyle="none", markersize=4.5, markerfacecolor=INK,
                          markeredgecolor=SURFACE, label="gauges only (no anchor)"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.125, 0.855), frameon=False, fontsize=6.2, ncol=6)
    title(fig, "What the interpolation and fusion methods predict on gauge-dry hours",
          "Held-out gauges, balanced spatial CV, 13,471 common hours; each anchor on its own evaluation mask. "
          "Gauge-only methods do not read a satellite and are drawn once, in ink.\nRain per year of scored hours: the grid "
          "follows FY4B's coverage and is weighted to summer, so compare methods, not with the full-record rates.", left=0.13)
    save(fig, "norain_fig11_fusion_dry")


def fig12_fusion_where(paired: list[dict], strata: list[dict], zeroing: list[dict]) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(10.0, 3.4), gridspec_kw={"width_ratios": [1.25, 1.0, 1.0, 1.0]})
    fig.subplots_adjust(left=0.135, right=0.99, top=0.72, bottom=0.2, wspace=0.42)
    methods = ["mgwr", "auto", "blend_mgwr", "blend_agrenv_mgwr"]
    ax = axes[0]
    clean_axis(ax, "x")
    quadrants = ["gauge dry, anchor dry", "gauge dry, anchor wet", "gauge wet, anchor dry", "gauge wet, anchor wet"]
    q_colors = ["#b7d4f6", "#256abf", "#f0c9a8", "#c35a1c"]
    anchor = "GPM"
    y = np.arange(len(methods))[::-1].astype(float)
    for k, method in enumerate(methods):
        rows = {r["level"]: r for r in pick(paired, scheme=SCHEME, anchor=anchor, method=method, stratifier="quadrant")}
        total = sum(rows[q]["sse_gap"] for q in quadrants)
        left_pos, left_neg = 0.0, 0.0
        for q, color in zip(quadrants, q_colors):
            gap = rows[q]["sse_gap"] / 1e3
            start = left_pos if gap >= 0 else left_neg
            ax.barh(y[k], gap, left=start, color=color, height=0.6, zorder=2, edgecolor=SURFACE, linewidth=0.8)
            if gap >= 0:
                left_pos += gap
            else:
                left_neg += gap
        ax.plot([total / 1e3], [y[k]], marker="|", markersize=11, color=INK, zorder=3)
    ax.axvline(0, color=AXIS, linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([METHOD_LABELS[m] for m in methods], fontsize=6.2)
    ax.set_xlabel("Squared-error gap to ADW (10³ mm², + = worse)")
    ax.set_title(f"Where the gap to ADW comes from · {anchor}", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in q_colors] + [Line2D([], [], marker="|", linestyle="none",
                                                                                  color=INK, markersize=8)]
    ax.legend(handles, quadrants + ["net gap"], frameon=False, fontsize=5.6, loc="upper center",
              bbox_to_anchor=(0.45, -0.2), ncol=3)

    for ax, stratifier, heading in ((axes[1], "distance", "Dry-hour mean by distance\nto the nearest training gauge"),
                                    (axes[2], "rain_proximity", "Dry-hour mean by distance\nto the gauge's own rain")):
        clean_axis(ax)
        levels = [r["level"] for r in pick(strata, anchor="GPM", method="adw", stratifier=stratifier, threshold=0.1)]
        x = np.arange(len(levels))
        for method, color, style in (("raw", MUTED, ":"), ("adw", INK, "--"), ("auto", ORDINAL_BLUES[1], "-"),
                                     ("blend_agrenv_mgwr", ORDINAL_BLUES[3], "-")):
            rows = [one(strata, anchor="GPM", method=method, stratifier=stratifier, level=lv, threshold=0.1) for lv in levels]
            ax.plot(x, [r["mean_dry"] for r in rows], color=color, linestyle=style, linewidth=1.5, marker="o",
                    markersize=3, label=METHOD_LABELS[method], zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels(levels, fontsize=5.8, rotation=30, ha="right")
        # Distance to rain spans more than a decade; distance to the training gauges does not.
        if stratifier == "rain_proximity":
            ax.set_yscale("log")
            for axis_formatter in (ax.yaxis.set_major_formatter, ax.yaxis.set_minor_formatter):
                axis_formatter(mpl.ticker.FuncFormatter(
                    lambda v, _: f"{v:g}" if v > 0 and round(v / 10 ** np.floor(np.log10(v))) in (1, 2, 5) else ""))
        ax.set_ylabel("mm/h on gauge-dry hours")
        ax.set_title(heading, loc="left")
    axes[2].legend(frameon=False, fontsize=5.6, loc="upper right")

    ax = axes[3]
    clean_axis(ax)
    for method, color, style in (("raw", MUTED, ":"), ("adw", INK, "--"), ("auto", ORDINAL_BLUES[1], "-"),
                                 ("blend_agrenv_mgwr", ORDINAL_BLUES[3], "-")):
        rows = sorted(pick(zeroing, anchor="GPM", method=method), key=lambda r: r["tau"])
        base = rows[0]["RMSE"]
        ax.plot([r["tau"] for r in rows], [100 * (base - r["RMSE"]) / base for r in rows], color=color, linestyle=style,
                linewidth=1.5, marker="o", markersize=3, label=METHOD_LABELS[method], zorder=2)
    ax.axhline(0, color=AXIS, linewidth=0.9)
    ax.set_xlabel("Values below τ set to 0 (mm/h)")
    ax.set_ylabel("Overall RMSE change (%, + = better)")
    ax.set_title("Zeroing small values · GPM\n(descriptive, not tuned in-fold)", loc="left")
    title(fig, "Where the methods' dry-hour error sits, and what zeroing drizzle would do",
          "Balanced spatial CV on the GPM anchor. Left: each gauge × anchor quadrant's contribution to the method's total "
          "squared error minus ADW's on the same cells.", left=0.135)
    save(fig, "norain_fig12_fusion_where")


# --------------------------------------------------------------------------------------------------
# Report tables


def report_tables(occurrence, exceedance, baseline, extent, network, hourly_summary, daily, daily_summary, lags,
                  zero_pattern, summary, paired, strata, fusion_days, zeroing, provenance, settings) -> None:
    overview = []
    for sample, products in SAMPLES.items():
        for product in products:
            for threshold in (0.1, 0.5):
                r = occ(occurrence, sample, product, "all", "all", threshold)
                overview.append({
                    "sample": sample, "product": product, "threshold": f"{threshold:g}",
                    "n_dry": f"{int(r['n_dry']):,}", "dry_share": pct(r["dry_share"]),
                    "POFD": pct(r["POFD"]), "POFD_ci": interval(r["POFD_lo"], r["POFD_hi"]),
                    "specificity": pct(r["specificity"]), "FAR": pct(r["FAR"], 0), "FAR_ci": interval(r["FAR_lo"], r["FAR_hi"], digits=0),
                    "NPV": pct(r["NPV"]), "NPV_ci": interval(r["NPV_lo"], r["NPV_hi"]),
                    "freq_bias": num(r["freq_bias"]), "freq_bias_ci": interval(r["freq_bias_lo"], r["freq_bias_hi"], 1, 2),
                    "PSS": num(r["PSS"]), "HSS": num(r["HSS"]), "HSS_ci": interval(r["HSS_lo"], r["HSS_hi"], 1, 2),
                    "POD": pct(r["POD"], 0),
                    "mean_dry": num(r["mean_dry"], 3), "spurious_mm": f"{r['spurious_mm_per_year']:.0f}",
                    "estimate_mm": f"{r['estimate_mm_per_year']:.0f}", "gauge_mm": f"{r['gauge_mm_per_year']:.0f}",
                    "dry_volume_share": pct(r["dry_volume_share"], 0), "zero_share": pct(r["zero_share_dry"], 0),
                })
    write_table("overview", overview)

    anatomy = []
    for stratifier in ("all", "season", "month", "region", "rain_proximity", "rain_side", "neighbours", "network",
                       "dewpoint_depression", "t2m"):
        rows = pick(occurrence, sample="gpm_gsmap_full", product="GPM", stratifier=stratifier, threshold=0.1)
        total_dry = sum(r["n_dry"] for r in rows)
        for r in rows:
            anatomy.append({"stratifier": stratifier, "level": r["level"], "dry_share": pct(r["dry_share"]),
                            "share_of_dry": pct(r["n_dry"] / total_dry) if total_dry else "–",
                            "gauge_mm": f"{r['gauge_mm_per_year']:.0f}"})
    hours = pick(occurrence, sample="gpm_gsmap_full", product="GPM", stratifier="hour", threshold=0.1)
    driest = max(hours, key=lambda r: r["dry_share"])
    wettest = min(hours, key=lambda r: r["dry_share"])
    anatomy.append({"stratifier": "hour_extremes", "level": "driest", "dry_share": pct(driest["dry_share"]),
                    "share_of_dry": driest["level"], "gauge_mm": ""})
    anatomy.append({"stratifier": "hour_extremes", "level": "wettest", "dry_share": pct(wettest["dry_share"]),
                    "share_of_dry": wettest["level"], "gauge_mm": ""})
    net = network[0]
    anatomy.append({"stratifier": "network_hours", "level": "network dry", "dry_share": pct(net["n_network_dry"] / net["n_hours_reporting"]),
                    "share_of_dry": f"{int(net['n_network_dry']):,}", "gauge_mm": f"{int(net['n_hours_reporting']):,}"})
    write_table("anatomy", anatomy)

    strata_rows = []
    for sample, products in SAMPLES.items():
        for product in products:
            for stratifier in ("season", "region", "rain_proximity", "rain_side", "neighbours", "network",
                               "dewpoint_depression", "t2m", "other_products"):
                for r05 in pick(occurrence, sample=sample, product=product, stratifier=stratifier, threshold=0.5):
                    r01 = occ(occurrence, sample, product, stratifier, r05["level"], 0.1)
                    strata_rows.append({
                        "sample": sample, "product": product, "stratifier": stratifier, "level": r05["level"],
                        "POFD01": pct(r01["POFD"]), "POFD01_ci": interval(r01["POFD_lo"], r01["POFD_hi"]),
                        "POFD05": pct(r05["POFD"]), "POFD05_ci": interval(r05["POFD_lo"], r05["POFD_hi"]),
                        "FAR01": pct(r01["FAR"], 0), "freq_bias01": num(r01["freq_bias"]),
                        "mean_dry": num(r05["mean_dry"], 3), "spurious_mm": f"{r05['spurious_mm_per_year']:.0f}",
                        "dry_volume_share": pct(r05["dry_volume_share"], 0),
                        "share_of_dry": pct(r05["n_dry"] / sum(x["n_dry"] for x in pick(
                            occurrence, sample=sample, product=product, stratifier=stratifier, threshold=0.5))),
                        "dry_share": pct(r05["dry_share"]),
                    })
    write_table("strata", strata_rows)

    hour_rows = []
    for sample, products in SAMPLES.items():
        for product in products:
            for stratifier in ("hour", "season_hour"):
                rows = pick(occurrence, sample=sample, product=product, stratifier=stratifier, threshold=0.5)
                groups = {}
                for r in rows:
                    season = r["level"].split("|")[0] if stratifier == "season_hour" else "all"
                    groups.setdefault(season, []).append(r)
                for season, members in groups.items():
                    peak = max(members, key=lambda r: r["POFD"])
                    low = min(members, key=lambda r: r["POFD"])
                    hour_rows.append({"sample": sample, "product": product, "season": season,
                                      "peak_hour": peak["level"].split("|")[-1], "peak_POFD05": pct(peak["POFD"]),
                                      "low_hour": low["level"].split("|")[-1], "low_POFD05": pct(low["POFD"]),
                                      "ratio": num(peak["POFD"] / low["POFD"], 1) if low["POFD"] > 0 else "–"})
    write_table("diurnal", hour_rows)

    station_rows = []
    for sample, products in SAMPLES.items():
        for product in products:
            for threshold in (0.1, 0.5):
                values = np.array([r["POFD"] for r in pick(occurrence, sample=sample, product=product, stratifier="station",
                                                            threshold=threshold)])
                spurious = np.array([r["spurious_mm_per_year"] for r in pick(occurrence, sample=sample, product=product,
                                                                              stratifier="station", threshold=threshold)])
                station_rows.append({"sample": sample, "product": product, "threshold": f"{threshold:g}",
                                     "POFD_p10": pct(np.quantile(values, 0.1)), "POFD_median": pct(np.median(values)),
                                     "POFD_p90": pct(np.quantile(values, 0.9)), "POFD_max": pct(values.max()),
                                     "spurious_median": f"{np.median(spurious):.0f}", "spurious_p90": f"{np.quantile(spurious, 0.9):.0f}",
                                     "spurious_max": f"{spurious.max():.0f}"})
    write_table("stations", station_rows)

    write_table("baseline", [{
        "sample": r["sample"], "level": r["level"], "n_pairs": f"{int(r['n_pairs']):,}", "median_km": num(r["median_km"], 1),
        "POFD": pct(r["POFD"]), "POFD_ci": interval(r["POFD_lo"], r["POFD_hi"]), "FAR": pct(r["FAR"], 0),
        "HSS": num(r["HSS"]), "POD": pct(r["POD"], 0), "mean_dry": num(r["mean_dry"], 3),
    } for r in baseline])

    write_table("extent", [{
        "sample": r["sample"], "product": r["product"], "n_hours": f"{int(r['n_hours']):,}", "mean_wet_share": pct(r["mean_wet_share"]),
        "any": pct(r["share_hours_gt_0"], 0), "gt5": pct(r["share_hours_gt_5"], 0), "gt25": pct(r["share_hours_gt_25"], 0),
        "gt50": pct(r["share_hours_gt_50"], 1),
    } for r in extent])

    write_table("spells", [{
        "sample": r["sample"], "source": r["source"], "threshold": f"{r['threshold']:g}", "n_spells": f"{int(r['n_spells']):,}",
        "dry_share": pct(r["dry_share"]), "mean_length": num(r["mean_length"], 1), "median_length": num(r["median_length"], 0),
        "p90_length": num(r["p90_length"], 0), "max_length": f"{int(r['max_length'])}",
        "share_in_ge24": pct(r["share_in_ge24"], 0), "share_in_ge72": pct(r["share_in_ge72"], 0),
    } for r in hourly_summary if r["level"] == "all"])

    day_rows = []
    for r in daily:
        if r["stratifier"] not in ("all", "season"):
            continue
        day_rows.append({"sample": r["sample"], "product": r["product"], "stratifier": r["stratifier"], "level": r["level"],
                         "threshold": f"{r['threshold']:g}", "n_days": f"{int(r['n_cells']):,}", "dry_share": pct(r["dry_share"]),
                         "POFD": pct(r["POFD"]), "POFD_ci": interval(r["POFD_lo"], r["POFD_hi"]), "NPV": pct(r["NPV"]),
                         "FAR": pct(r["FAR"], 0), "POD": pct(r["POD"], 0), "freq_bias": num(r["freq_bias"]),
                         "freq_bias_ci": interval(r["freq_bias_lo"], r["freq_bias_hi"], 1, 2), "HSS": num(r["HSS"]),
                         "spurious_mm": f"{r['spurious_mm_per_year'] / 24:.0f}"})
    write_table("days", day_rows)

    cdd = []
    for r in daily_summary:
        if r["level"] != "all":
            continue
        stations = [x for x in pick(daily_summary, sample=r["sample"], pair=r["pair"], source=r["source"],
                                    threshold=r["threshold"]) if x["level"] != "all"]
        cdd.append({"sample": r["sample"], "pair": r["pair"], "source": r["source"], "threshold": f"{r['threshold']:g}",
                    "mean_length": num(r["mean_length"], 1), "max_length": f"{int(r['max_length'])}",
                    "median_station_max": f"{np.median([x['max_length'] for x in stations]):.0f}",
                    "dry_share": pct(r["dry_share"]), "n_spells": f"{int(r['n_spells']):,}"})
    write_table("dry_spells_daily", cdd)

    calibration = [{"kind": "lag", "lag": f"{int(r['lag'])}", "r": num(r["r"], 3)} for r in lags]
    z = zero_pattern[0]
    calibration.append({"kind": "pattern", **{k: pct(v, 1) if k != "n" and "ratio" not in k else (
        f"{int(v):,}" if k == "n" else num(v, 2)) for k, v in z.items()}})
    write_table("calibration", calibration)

    fusion_rows = []
    for r in summary:
        if r["threshold"] != 0.1:
            continue
        r05 = one(summary, scheme=r["scheme"], anchor=r["anchor"], method=r["method"], stratifier="all", threshold=0.5)
        dry = pick(paired, scheme=r["scheme"], anchor=r["anchor"], method=r["method"], stratifier="gauge_state", level="gauge dry")
        quad = {x["level"]: x for x in pick(paired, scheme=r["scheme"], anchor=r["anchor"], method=r["method"], stratifier="quadrant")}
        fusion_rows.append({
            "scheme": r["scheme"], "anchor": r["anchor"], "method": r["method"],
            "n_dry": f"{int(r['n_dry']):,}", "mean_dry": num(r["mean_dry"], 3), "mean_dry_ci": interval(r["mean_dry_lo"], r["mean_dry_hi"], 1, 3),
            "spurious_mm": f"{r['spurious_mm_per_year']:.0f}", "gauge_mm": f"{r['gauge_mm_per_year']:.0f}",
            "POFD01": pct(r["POFD"]), "POFD05": pct(r05["POFD"]), "zero_share": pct(r["zero_share_dry"], 0),
            "RMSE_dry": num(r["RMSE_dry"], 3), "dry_sse_share": pct(r["dry_sse_share"], 0), "FAR01": pct(r["FAR"], 0),
            "freq_bias01": num(r["freq_bias"]), "HSS01": num(r["HSS"]),
            "dry_vs_adw": pct(dry[0]["improvement"]) if dry else "",
            "dry_vs_adw_ci": interval(dry[0]["improvement_lo"], dry[0]["improvement_hi"]) if dry else "",
            "gap_dry_anchor_wet": pct(quad["gauge dry, anchor wet"]["sse_gap_share"], 0) if quad else "",
            "gap_dry_anchor_dry": pct(quad["gauge dry, anchor dry"]["sse_gap_share"], 0) if quad else "",
            "gap_dry_total": pct(quad["gauge dry, anchor wet"]["sse_gap_share"] + quad["gauge dry, anchor dry"]["sse_gap_share"], 0) if quad else "",
            "gap_total": num(sum(x["sse_gap"] for x in quad.values()) / 1e3, 1) if quad else "",
        })
    write_table("fusion", fusion_rows)

    fusion_strata = []
    for r in strata:
        if r["threshold"] != 0.1:
            continue
        fusion_strata.append({"anchor": r["anchor"], "method": r["method"], "stratifier": r["stratifier"], "level": r["level"],
                              "mean_dry": num(r["mean_dry"], 3), "POFD01": pct(r["POFD"]), "RMSE_dry": num(r["RMSE_dry"], 3),
                              "zero_share": pct(r["zero_share_dry"], 0), "share_of_dry": "",
                              "n_dry": f"{int(r['n_dry']):,}"})
    write_table("fusion_strata", fusion_strata)

    write_table("fusion_zero", [{
        "anchor": r["anchor"], "method": r["method"], "tau": f"{r['tau']:g}", "RMSE": num(r["RMSE"], 4),
        "RMSE_dry": num(r["RMSE_dry"], 3), "mean_dry": num(r["mean_dry"], 3), "POD": pct(r["POD"], 0),
        "POFD": pct(r["POFD"]), "CSI": num(r["CSI"], 3), "volume_ratio": num(r["volume_ratio"], 2),
        "RMSE_change": pct((min(x["RMSE"] for x in pick(zeroing, anchor=r["anchor"], method=r["method"], tau=0.0)) - r["RMSE"]) /
                           min(x["RMSE"] for x in pick(zeroing, anchor=r["anchor"], method=r["method"], tau=0.0)), 1),
    } for r in zeroing])

    write_table("fusion_days", [{
        "anchor": r["anchor"], "method": r["method"], "stratifier": r["stratifier"], "level": r["level"],
        "threshold": f"{r['threshold']:g}", "POFD": pct(r["POFD"]), "NPV": pct(r["NPV"]), "FAR": pct(r["FAR"], 0),
        "POD": pct(r["POD"], 0), "freq_bias": num(r["freq_bias"]), "HSS": num(r["HSS"]), "dry_share": pct(r["dry_share"]),
        "n_days": f"{int(r['n_cells']):,}",
    } for r in fusion_days])

    write_table("provenance", [{"key": r["key"], "value": r["value"]} for r in provenance] +
                [{"key": f"setting_{r['setting']}", "value": str(r["value"])} for r in settings])


def silent_tables(silent, sensitivity, occurrence, paired) -> None:
    flagged = sorted([r for r in silent if r["flagged"]], key=lambda r: -r["neighbour_total"])
    write_table("silent_spells", [{
        "station_id": r["station_id"], "first_date": r["first_date"], "last_date": r["last_date"], "season": r["season"],
        "n_days": f"{int(r['n_days'])}", "gauge_total": f"{r['gauge_total']:.1f}", "neighbour_total": f"{r['neighbour_total']:.0f}",
        "neighbour_wet_days": f"{int(r['neighbour_wet_days'])}", "GPM_total": f"{r['GPM_total']:.0f}",
        "GSMaP_total": f"{r['GSMaP_total']:.0f}",
    } for r in flagged])
    seasons = {}
    for r in flagged:
        seasons[r["season"]] = seasons.get(r["season"], 0) + int(r["n_days"])
    total_days = sum(seasons.values())
    write_table("silent_seasons", [{"season": s, "station_days": f"{seasons.get(s, 0):,}",
                                    "share": pct(seasons.get(s, 0) / total_days, 0)} for s in SEASONS])
    write_table("silent_sensitivity", [{
        "min_wet_days": f"{int(r['min_wet_days'])}", "n_spells": f"{int(r['n_spells'])}", "n_stations": f"{int(r['n_stations'])}",
        "station_days": f"{int(r['station_days']):,}", "dry_hours": f"{int(r['dry_hours_in_spells']):,}",
        "dry_hour_share": pct(r["dry_hour_share"], 1), "p_dry_given_wet": pct(r["p_dry_given_wet"], 1),
    } for r in sensitivity])
    screened = []
    for sample, products in SAMPLES.items():
        for product in products:
            for threshold in (0.1, 0.5):
                base = occ(occurrence, sample, product, "all", "all", threshold)
                ok = occ(occurrence, sample, product, "gauge_screen", "gauge working", threshold)
                bad = occ(occurrence, sample, product, "gauge_screen", "silent-gauge spell", threshold)
                screened.append({
                    "sample": sample, "product": product, "threshold": f"{threshold:g}", "season": "all",
                    "POFD_all": pct(base["POFD"]), "POFD_working": pct(ok["POFD"]),
                    "POFD_working_ci": interval(ok["POFD_lo"], ok["POFD_hi"]), "POFD_silent": pct(bad["POFD"]),
                    "POFD_silent_ci": interval(bad["POFD_lo"], bad["POFD_hi"]),
                    "FAR_all": pct(base["FAR"], 0), "FAR_working": pct(ok["FAR"], 0),
                    "freq_bias_all": num(base["freq_bias"]), "freq_bias_working": num(ok["freq_bias"]),
                    "silent_dry_share": pct(bad["n_dry"] / base["n_dry"], 1),
                    "silent_false_alarm_share": pct(bad["n_dry"] * bad["POFD"] / (base["n_dry"] * base["POFD"]), 0),
                    "spurious_all": f"{base['spurious_mm_per_year']:.0f}", "spurious_working": f"{ok['spurious_mm_per_year']:.0f}",
                })
                for season in SEASONS:
                    b = occ(occurrence, sample, product, "season", season, threshold)
                    w = occ(occurrence, sample, product, "season_screened", season, threshold)
                    screened.append({"sample": sample, "product": product, "threshold": f"{threshold:g}", "season": season,
                                     "POFD_all": pct(b["POFD"]), "POFD_working": pct(w["POFD"]),
                                     "FAR_all": pct(b["FAR"], 0), "FAR_working": pct(w["FAR"], 0),
                                     "freq_bias_all": num(b["freq_bias"]), "freq_bias_working": num(w["freq_bias"])})
    write_table("screened", screened)
    shares = []
    for r in paired:
        if r["stratifier"] != "screened_state" or r["scheme"] != SCHEME:
            continue
        shares.append({"anchor": r["anchor"], "method": r["method"], "level": r["level"], "n": f"{int(r['n']):,}",
                       "gap_share": pct(r["sse_gap_share"], 0), "improvement": pct(r["improvement"]),
                       "improvement_ci": interval(r["improvement_lo"], r["improvement_hi"]),
                       "mean_method": num(r["mean_method"], 3), "mean_ref": num(r["mean_ref"], 3)})
    write_table("fusion_screened", shares)


def main() -> None:
    occurrence = read_rows(RUN_DIR / "dry_occurrence.csv")
    exceedance = read_rows(RUN_DIR / "dry_exceedance.csv")
    baseline = read_rows(RUN_DIR / "neighbour_gauge_baseline.csv")
    extent = read_rows(RUN_DIR / "phantom_extent.csv")
    network = read_rows(RUN_DIR / "network_hours.csv")
    hourly_summary = read_rows(RUN_DIR / "dry_spells_hourly_summary.csv")
    hourly_hist = read_rows(RUN_DIR / "dry_spells_hourly_histogram.csv")
    daily = read_rows(RUN_DIR / "dry_days.csv")
    daily_summary = read_rows(RUN_DIR / "dry_spells_daily_summary.csv")
    lags = read_rows(RUN_DIR / "gpm_calibration_lags.csv")
    zero_pattern = read_rows(RUN_DIR / "gpm_calibration_zero_pattern.csv")
    stations = read_rows(RUN_DIR / "stations.csv")
    settings = read_rows(RUN_DIR / "run_settings.csv")
    silent = read_rows(RUN_DIR / "silent_gauge_spells.csv")
    sensitivity = read_rows(RUN_DIR / "silent_gauge_sensitivity.csv")
    summary = read_rows(FUSION_DIR / "fusion_dry_summary.csv")
    paired = read_rows(FUSION_DIR / "fusion_dry_paired.csv")
    strata = read_rows(FUSION_DIR / "fusion_dry_strata.csv")
    fusion_days = read_rows(FUSION_DIR / "fusion_dry_days.csv")
    zeroing = read_rows(FUSION_DIR / "fusion_zero_threshold.csv")
    provenance = read_rows(FUSION_DIR / "source_provenance.csv")
    for row in settings:
        row["value"] = str(row["value"])

    fig1_anatomy(occurrence)
    fig2_exceedance(exceedance, baseline)
    fig3_contingency(occurrence)
    fig4_diurnal(occurrence)
    fig5_maps(occurrence, stations)
    fig6_context(occurrence)
    fig7_floor_and_air(occurrence, baseline)
    fig8_agreement(occurrence, extent)
    fig9_spells(hourly_hist, daily, daily_summary)
    fig10_calibration(occurrence, exceedance, lags)
    fig11_fusion(summary, paired)
    fig12_fusion_where(paired, strata, zeroing)
    fig13_silent(silent, occurrence)

    report_tables(occurrence, exceedance, baseline, extent, network, hourly_summary, daily, daily_summary, lags,
                  zero_pattern, summary, paired, strata, fusion_days, zeroing, provenance, settings)
    silent_tables(silent, sensitivity, occurrence, paired)
    figures = REPORT_DIR / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    for png in sorted(FIG_DIR.glob("norain_*.png")):
        shutil.copy2(png, figures / png.name)
    print(f"Copied figures to {figures}")


if __name__ == "__main__":
    main()
