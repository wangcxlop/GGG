"""Draw the satellite temporal-evaluation figures.

Reads the CSVs written by `scripts/run_satellite_temporal_evaluation.jl` (and the landform diagnostics
of `scripts/prepare_station_landform.jl`) and writes PNG (300 dpi) and PDF figures to
`output/satellite_temporal_evaluation/figures/`:

    py -3.13 scripts/plot_satellite_temporal_evaluation.py

Two samples run through every figure: `all_products` (the hours FY4B covers, so all three products on
identical cells) and `gpm_gsmap_full` (GPM and GSMaP over the whole 2022-2024 record). Regions with fewer
than 20 gauges are drawn hollow or left out of regional contrasts and are never interpreted.

Colours follow the validated palette of `plot_heavy_rain_events.py`: the first three categorical slots
for the products (validated all-pairs), gauge in primary ink, a single-hue blue ramp for magnitudes and
ordered classes, and a blue/red diverging ramp with a grey centre for signed errors.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "output" / "satellite_temporal_evaluation"
FIG_DIR = RUN_DIR / "figures"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
PRODUCT_COLORS = {"FY4B": "#2a78d6", "GPM": "#eb6834", "GSMaP": "#1baf7a"}
SAMPLE_PRODUCTS = {"all_products": ["FY4B", "GPM", "GSMaP"], "gpm_gsmap_full": ["GPM", "GSMaP"]}
SAMPLE_TITLES = {
    "all_products": "All three products · hours FY4B covers",
    "gpm_gsmap_full": "GPM and GSMaP · full 2022–2024 record",
}
SEASONS = ["all", "MAM", "JJA", "SON", "DJF"]
CLASSES = ["light", "moderate", "heavy", "rainstorm", "severe_rainstorm"]
CLASS_LABELS = ["Light\n0.1–2", "Moderate\n2–4", "Heavy\n4–8", "Rainstorm\n8–20", "Severe\n≥20"]
REGIONS = ["relief_lt500", "relief_500_1000", "relief_ge1000"]
REGION_LABELS = {"all": "All gauges", "relief_lt500": "Relief < 500 m", "relief_500_1000": "Relief 500–1000 m",
                 "relief_ge1000": "Relief ≥ 1000 m"}
REGION_MARKERS = {"relief_lt500": "o", "relief_500_1000": "s", "relief_ge1000": "^"}
# Ordinal blue ramp (validated --ordinal): lightest step still clears 2:1 on the surface.
ORDINAL_BLUES = ["#86b6ef", "#5598e7", "#256abf", "#104281"]
SEQUENTIAL_BLUES = ["#f0efec", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
CAVEAT = ("Gauges are 0.5 mm tipping buckets: gauge 'light' hours are 0.5–1.5 mm and satellite 0.1–0.5 mm "
          "cannot be verified.\nHollow markers: fewer than 100 hours or 30 events. Winter gauge timing may include "
          "delayed snowmelt tips (the DJF gauge diurnal cycle peaks at midday).")

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


TEXT_COLUMNS = {"station_id"}   # numeric-looking identifiers stay text


def read_rows(path: Path) -> list[dict]:
    """CSV rows with numeric-looking fields converted to float and true/false to bool."""
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            parsed = {}
            for key, value in row.items():
                if key in TEXT_COLUMNS:
                    parsed[key] = value
                    continue
                if value in ("true", "false"):
                    parsed[key] = value == "true"
                    continue
                try:
                    parsed[key] = float(value)
                except ValueError:
                    parsed[key] = value
            rows.append(parsed)
    return rows


def pick(rows: list[dict], **conditions) -> list[dict]:
    return [row for row in rows if all(row.get(k) == v for k, v in conditions.items())]


def one(rows: list[dict], **conditions) -> dict | None:
    found = pick(rows, **conditions)
    return found[0] if found else None


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


def product_legend(fig: plt.Figure, products: list[str], extra: list | None = None, **kwargs) -> None:
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=5.5, markerfacecolor=PRODUCT_COLORS[p],
                      markeredgecolor=SURFACE, label=p) for p in products]
    handles += extra or []
    fig.legend(handles=handles, frameon=False, fontsize=6.8, ncol=len(handles), **kwargs)


def dots(ax, x, values, color, lo=None, hi=None, hollow=None, size=22):
    """One product's values as dots with optional interval whiskers; `hollow` marks low-sample points."""
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    hollow = np.zeros(len(values), dtype=bool) if hollow is None else np.asarray(hollow, dtype=bool)
    if lo is not None:
        ax.vlines(x, lo, hi, color=color, linewidth=1.0, alpha=0.7, zorder=2)
    ax.scatter(x[~hollow], values[~hollow], s=size, color=color, edgecolors=SURFACE, linewidths=0.8, zorder=3)
    ax.scatter(x[hollow], values[hollow], s=size, facecolors=SURFACE, edgecolors=color, linewidths=1.0, zorder=3)


# --------------------------------------------------------------------------------------------------
# Figure 1: landform regions


def fig1_landform() -> None:
    landform_dir = RUN_DIR / "landform"
    blocks = read_rows(landform_dir / "relief_blocks.csv")
    change = read_rows(landform_dir / "window_change_point.csv")
    stations = read_rows(ROOT / "data" / "processed" / "covariates" / "station_landform.csv")
    selected = next(row for row in change if row["selected"])

    fig = plt.figure(figsize=(7.2, 5.0))
    ax_map = fig.add_axes([0.05, 0.22, 0.44, 0.64])
    ax_relief = fig.add_axes([0.6, 0.6, 0.37, 0.26])
    ax_gain = fig.add_axes([0.6, 0.22, 0.37, 0.26])

    cmap = ListedColormap(ORDINAL_BLUES[:3])
    norm = BoundaryNorm([0, 500, 1000, 5000], cmap.N)
    lon = np.array([b["lon"] for b in blocks])
    lat = np.array([b["lat"] for b in blocks])
    relief = np.array([b["relief_m"] for b in blocks])
    ax_map.scatter(lon, lat, c=relief, cmap=cmap, norm=norm, marker="s", s=5.5, linewidths=0, zorder=1, rasterized=True)
    for region in REGIONS:
        members = [s for s in stations if s["relief_region"] == region]
        insufficient = members and members[0]["insufficient_gauges"]
        ax_map.scatter([s["lon"] for s in members], [s["lat"] for s in members], marker=REGION_MARKERS[region], s=11,
                       facecolors=SURFACE if insufficient else INK, edgecolors=INK, linewidths=0.6, zorder=3)
    mid_lat = (lat.min() + lat.max()) / 2
    ax_map.set_aspect(1 / np.cos(np.radians(mid_lat)), adjustable="box")
    ax_map.set_xticks([110, 110.5, 111])
    ax_map.set_xticklabels(["110°E", "110.5°E", "111°E"])
    ax_map.set_yticks([31.5, 32, 32.5, 33])
    ax_map.set_yticklabels(["31.5°N", "32°N", "32.5°N", "33°N"])
    ax_map.tick_params(length=2, pad=1.5, labelsize=6)
    counts = {r: sum(s["relief_region"] == r for s in stations) for r in REGIONS}
    handles = [Patch(facecolor=ORDINAL_BLUES[k], edgecolor="none",
                     label=f"{REGION_LABELS[r]} · {counts[r]} gauges" + (" (insufficient)" if counts[r] < 20 else ""))
               for k, r in enumerate(REGIONS)]
    handles += [Line2D([], [], linestyle="none", marker=REGION_MARKERS[r], markersize=4, markerfacecolor=INK if counts[r] >= 20 else SURFACE,
                       markeredgecolor=INK, label=f"Gauge, {REGION_LABELS[r].lower()}") for r in REGIONS]
    ax_map.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2, frameon=False, fontsize=5.8)
    ax_map.set_title(f"Local relief in {selected['side_km']:g} km windows", loc="left", pad=4)

    side = [row["side_km"] for row in change]
    for ax, key, label in ((ax_relief, "mean_relief_m", "Area-wide mean relief (m)"),
                           (ax_gain, "change_point_gain", "Change-point gain S − Sᵢ")):
        clean_axis(ax)
        values = [row[key] for row in change]
        ax.plot(side, values, color=PRODUCT_COLORS["FY4B"], linewidth=1.6, zorder=2)
        ax.axvline(selected["side_km"], color=AXIS, linewidth=0.9, zorder=1)
        ax.scatter([selected["side_km"]], [selected[key]], s=26, color=INK, edgecolors=SURFACE, linewidths=0.8, zorder=3)
        ax.set_title(label, loc="left", pad=3)
    ax_gain.set_xlabel("Window side (km)")
    ax_gain.annotate(f"{selected['side_km']:g} km", (selected["side_km"], selected["change_point_gain"]),
                     xytext=(6, -10), textcoords="offset points", fontsize=6, color=INK)
    fig.suptitle("Landform regions from Copernicus GLO-30 relief (no gauge sits on a plain)", x=0.05, ha="left",
                 fontsize=9, y=0.975)
    fig.text(0.05, 0.915, "Window fixed by the mean change-point of ln(relief / area) before classifying any gauge. "
             "Breaks: hills < 200 m, small-relief mountains < 500 m, medium < 1000 m.", ha="left", fontsize=6.3,
             color=INK_SECONDARY)
    save(fig, "fig1_landform_regions")


# --------------------------------------------------------------------------------------------------
# Figure 2: intensity classes

INTENSITY_PANELS = [
    ("RB_pct", "Relative bias (%)", 0.0, True),
    ("POD_rain", "Detected as rain (sat ≥ 0.1)", 1.0, False),
    ("class_hit", "Same intensity class", 1.0, True),
    ("under_class", "Satellite in a lower class", 0.0, False),
]


def fig2_intensity(metrics: list[dict], sample: str) -> None:
    products = SAMPLE_PRODUCTS[sample]
    fig, axes = plt.subplots(len(INTENSITY_PANELS), len(SEASONS), figsize=(10.0, 7.6), sharey="row")
    fig.subplots_adjust(left=0.07, right=0.975, top=0.84, bottom=0.11, wspace=0.08, hspace=0.42)
    x = np.arange(len(CLASSES))
    offsets = np.linspace(-0.22, 0.22, len(products))
    for j, season in enumerate(SEASONS):
        for i, (key, title, ideal, with_ci) in enumerate(INTENSITY_PANELS):
            ax = axes[i, j]
            clean_axis(ax)
            ax.axhline(ideal, color=AXIS, linewidth=0.9, zorder=1)
            for offset, product in zip(offsets, products):
                rows = [one(metrics, sample=sample, product=product, season=season, region="all", gauge_class=c) for c in CLASSES]
                values = [r[key] if r and r["n"] > 0 else np.nan for r in rows]
                lo = [r[f"{key}_lo"] if r and with_ci and r["n"] > 0 else np.nan for r in rows] if with_ci else None
                hi = [r[f"{key}_hi"] if r and with_ci and r["n"] > 0 else np.nan for r in rows] if with_ci else None
                hollow = [bool(r["low_sample"]) if r else True for r in rows]
                dots(ax, x + offset, values, PRODUCT_COLORS[product], lo, hi, hollow)
            ax.set_xticks(x)
            ax.set_xlim(-0.6, len(CLASSES) - 0.4)
            ax.set_xticklabels(CLASS_LABELS if i == len(INTENSITY_PANELS) - 1 else [], fontsize=5.6)
            if i == 0:
                ax.set_title("All seasons" if season == "all" else season, pad=4)
            if j == 0:
                ax.set_ylabel(title, fontsize=6.5)
    axes[-1, len(SEASONS) // 2].set_xlabel("Gauge hourly intensity class (mm/h)")
    product_legend(fig, products, loc="upper right", bbox_to_anchor=(0.99, 0.975))
    fig.suptitle(f"Hourly performance by gauge intensity class — {SAMPLE_TITLES[sample]}", x=0.07, ha="left", fontsize=9, y=0.975)
    fig.text(0.07, 0.935, "Gauge-wet hours only (no-rain hours excluded), all gauges. Whiskers: 95% day-block bootstrap. "
             "Grey line: perfect score.", ha="left", fontsize=6.3, color=INK_SECONDARY)
    fig.text(0.07, 0.918, CAVEAT, ha="left", va="top", fontsize=6.0, color=MUTED, linespacing=1.4)
    save(fig, f"fig2_intensity_classes_{sample}")


# --------------------------------------------------------------------------------------------------
# Figure 3: confusion matrices

SAT_CLASSES = ["no_rain"] + CLASSES
SAT_LABELS = ["No rain", "Light", "Moderate", "Heavy", "Rainstorm", "Severe"]


def fig3_confusion(confusion: list[dict]) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(8.4, 5.6))
    fig.subplots_adjust(left=0.17, right=0.9, top=0.86, bottom=0.1, wspace=0.12, hspace=0.42)
    cmap = ListedColormap(SEQUENTIAL_BLUES)
    norm = BoundaryNorm([0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0], cmap.N)
    for i, sample in enumerate(SAMPLE_PRODUCTS):
        for j in range(3):
            ax = axes[i, j]
            products = SAMPLE_PRODUCTS[sample]
            if j >= len(products):
                ax.axis("off")
                continue
            product = products[j]
            grid = np.array([[one(confusion, sample=sample, product=product, season="all", region="all",
                                  gauge_class=g, sat_class=s)["row_fraction"] for s in SAT_CLASSES] for g in CLASSES])
            ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto")
            for r in range(grid.shape[0]):
                for c in range(grid.shape[1]):
                    value = grid[r, c]
                    if np.isfinite(value) and value >= 0.005:
                        ax.text(c, r, f"{100 * value:.0f}", ha="center", va="center", fontsize=5.6,
                                color=SURFACE if value >= 0.3 else INK)
            for r in range(grid.shape[0]):
                ax.add_patch(plt.Rectangle((r + 1 - 0.5, r - 0.5), 1, 1, fill=False, edgecolor=INK, linewidth=0.7))
            ax.set_xticks(range(len(SAT_CLASSES)))
            ax.set_xticklabels(SAT_LABELS, rotation=40, ha="right", fontsize=5.8)
            ax.set_yticks(range(len(CLASSES)))
            ax.set_yticklabels([label.split("\n")[0] for label in CLASS_LABELS] if j == 0 else [], fontsize=5.8)
            ax.tick_params(length=0, pad=2)
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.set_title(f"{product}", loc="left", pad=3, color=INK)
        axes[i, 0].text(-0.36, 0.5, SAMPLE_TITLES[sample].replace(" · ", "\n"), transform=axes[i, 0].transAxes,
                        rotation=90, ha="center", va="center", fontsize=6.5)
    cax = fig.add_axes([0.92, 0.25, 0.012, 0.5])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[0, 0.1, 0.3, 0.5, 1.0])
    bar.ax.set_yticklabels(["0", "10", "30", "50", "100"])
    bar.outline.set_edgecolor(AXIS)
    bar.ax.tick_params(length=2, labelsize=6)
    bar.set_label("Share of the gauge class (%)", fontsize=6.5, color=INK_SECONDARY)
    fig.suptitle("Where each gauge intensity class lands in the satellite's classes", x=0.17, ha="left", fontsize=9, y=0.975)
    fig.text(0.17, 0.935, "Rows: gauge class (gauge-wet hours, all seasons and gauges); columns: satellite class; outlined "
             "cells: same class.\nCells right of the outline are overestimates, left of it underestimates.",
             ha="left", va="top", fontsize=6.3, color=INK_SECONDARY, linespacing=1.4)
    save(fig, "fig3_intensity_confusion")


# --------------------------------------------------------------------------------------------------
# Figure 4: false alarms

FALSE_ALARM_STACK = [
    ("dry_share_below_gauge_resolution", "0.1–0.5 mm (below gauge resolution)"),
    ("dry_share_light_resolved", "0.5–2 mm"),
    ("dry_share_moderate", "2–4 mm"),
    ("heavy_plus", "≥ 4 mm"),
]


def fig4_false_alarms(false_alarms: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.6), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.76, bottom=0.2, wspace=0.08)
    for ax, sample in zip(axes, SAMPLE_PRODUCTS):
        clean_axis(ax)
        products = SAMPLE_PRODUCTS[sample]
        width = 0.8 / len(products)
        for k, product in enumerate(products):
            for s, season in enumerate(SEASONS):
                row = one(false_alarms, sample=sample, product=product, season=season, region="all")
                if row is None or not np.isfinite(row["false_alarm_rate"]):
                    continue
                x = s - 0.4 + width * (k + 0.5)
                bottom = 0.0
                for m, (key, _) in enumerate(FALSE_ALARM_STACK):
                    value = (row["dry_share_heavy"] + row["dry_share_rainstorm"] + row["dry_share_severe_rainstorm"]
                             if key == "heavy_plus" else row[key]) * 100
                    ax.bar(x, value, width=width * 0.86, bottom=bottom, color=ORDINAL_BLUES[m], edgecolor=SURFACE,
                           linewidth=0.6, zorder=2)
                    bottom += value
                if s == 0:
                    ax.text(x, bottom, product, rotation=90, ha="center", va="bottom", fontsize=5.4, color=INK_SECONDARY)
        ax.set_xticks(range(len(SEASONS)))
        ax.set_xticklabels(["All" if s == "all" else s for s in SEASONS])
        ax.set_title(SAMPLE_TITLES[sample], loc="left", pad=4)
    axes[0].set_ylabel("Satellite wet on gauge-dry hours (%)")
    handles = [Patch(facecolor=ORDINAL_BLUES[m], edgecolor="none", label=label) for m, (_, label) in enumerate(FALSE_ALARM_STACK)]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=6.5, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("False alarms: how often, and how hard, the satellite rains when the gauge is dry", x=0.08, ha="left",
                 fontsize=9, y=0.975)
    fig.text(0.08, 0.9, "Bars within a season, left to right: " + " / ".join(SAMPLE_PRODUCTS["all_products"]) +
             " (left panel), GPM / GSMaP (right). Stack: the satellite's hourly amount on those hours.",
             ha="left", fontsize=6.3, color=INK_SECONDARY)
    fig.text(0.08, 0.85, "The lightest segment is below the gauges' 0.5 mm resolution and may be real rain the bucket "
             "had not yet tipped.", ha="left", fontsize=6.0, color=MUTED)
    save(fig, "fig4_false_alarms")


# --------------------------------------------------------------------------------------------------
# Figure 5: station events

EVENT_PANELS = [
    ("POD_event", "Event detected (POD)", 1.0, True),
    ("peak_within_1h", "Peak hour within ±1 h", 1.0, True),
    ("median_centroid_error_h", "Rain-centre timing error (h, + = late)", 0.0, False),
    ("lag_balance", "Best lag: late minus early share", 0.0, False),
    ("pooled_volume_rel_bias", "Event volume bias", 0.0, True),
]


def fig5_events(summary: list[dict]) -> None:
    classes = ["all"] + CLASSES
    labels = ["All"] + [label.split("\n")[0] for label in CLASS_LABELS]
    for row in summary:
        row["lag_balance"] = row["lag_late_share"] - row["lag_early_share"]
    fig, axes = plt.subplots(2, len(EVENT_PANELS), figsize=(12.0, 5.2))
    fig.subplots_adjust(left=0.08, right=0.99, top=0.8, bottom=0.13, wspace=0.3, hspace=0.45)
    x = np.arange(len(classes))
    for i, sample in enumerate(SAMPLE_PRODUCTS):
        products = SAMPLE_PRODUCTS[sample]
        offsets = np.linspace(-0.22, 0.22, len(products))
        for j, (key, title, ideal, with_ci) in enumerate(EVENT_PANELS):
            ax = axes[i, j]
            clean_axis(ax)
            ax.axhline(ideal, color=AXIS, linewidth=0.9, zorder=1)
            for offset, product in zip(offsets, products):
                rows = [one(summary, min_dry_gap_h=3.0, sample=sample, product=product, season="all", region="all",
                            event_class=c) for c in classes]
                values = [r[key] for r in rows]
                if with_ci:
                    lo = [r[f"{key}_lo"] for r in rows]
                    hi = [r[f"{key}_hi"] for r in rows]
                else:
                    lo = hi = None
                dots(ax, x + offset, values, PRODUCT_COLORS[product], lo, hi, [r["low_sample"] for r in rows])
            ax.set_xticks(x)
            ax.set_xticklabels(labels if i == 1 else [], rotation=35, ha="right", fontsize=5.8)
            if i == 0:
                ax.set_title(title, loc="left", pad=4)
        axes[i, 0].text(-0.45, 0.5, SAMPLE_TITLES[sample].replace(" · ", "\n"), transform=axes[i, 0].transAxes,
                        rotation=90, ha="center", va="center", fontsize=6.5)
    fig.text(0.55, 0.035, "Event class = the gauge's peak hourly intensity", ha="center", fontsize=6.5, color=INK_SECONDARY)
    product_legend(fig, ["FY4B", "GPM", "GSMaP"], loc="upper right", bbox_to_anchor=(0.99, 0.975))
    fig.suptitle("Station rain events: detection, timing and volume", x=0.08, ha="left", fontsize=9, y=0.975)
    fig.text(0.08, 0.93, "Events end after 3 dry hours; scored over the event ± 3 h. Whiskers: 95% day-block bootstrap. "
             "Timing panels over detected events: median difference of rain-weighted mean hours, and the balance of best "
             "lags within ± 3 h.", ha="left", fontsize=6.3, color=INK_SECONDARY)
    fig.text(0.08, 0.895, "Peak hours are hourly and tipping-bucket ties are common, so the peak-hour median is 0 h almost "
             "everywhere; the rain-centre error resolves sub-hour shifts. Hollow: fewer than 30 events.",
             ha="left", fontsize=6.0, color=MUTED)
    save(fig, "fig5_event_timing")


# --------------------------------------------------------------------------------------------------
# Figure 6: diurnal cycles

def fig6_diurnal(diurnal: list[dict], summary: list[dict], sample: str, by: str) -> None:
    """Diurnal amount and wet frequency by season (region all) or by region (season all)."""
    products = SAMPLE_PRODUCTS[sample]
    strata = SEASONS if by == "season" else ["all", "relief_lt500", "relief_500_1000"]
    height = 1.55 * len(strata) + 1.6
    fig, axes = plt.subplots(len(strata), 2, figsize=(6.8, height), sharex=True)
    fig.subplots_adjust(left=0.13, right=0.97, top=1 - 1.25 / height, bottom=0.55 / height, hspace=0.38, wspace=0.22)
    for i, stratum in enumerate(strata):
        season, region = (stratum, "all") if by == "season" else ("all", stratum)
        for j, (key, title, scale) in enumerate((("mean_amount", "Mean amount (mm/h)", 1.0),
                                                 ("wet_freq", "Wet-hour frequency (%)", 100.0))):
            ax = axes[i, j]
            clean_axis(ax)
            for source, color, width in [("Gauge", INK, 1.8)] + [(p, PRODUCT_COLORS[p], 1.5) for p in products]:
                rows = sorted(pick(diurnal, sample=sample, source=source, season=season, region=region), key=lambda r: r["hour"])
                hours = [r["hour"] for r in rows if r["n"] >= 100]
                values = [r[key] * scale for r in rows if r["n"] >= 100]
                ax.plot(hours, values, color=color, linewidth=width, zorder=3 if source == "Gauge" else 2)
            if i == 0:
                ax.text(0, 1.2, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, color=INK)
            if j == 0:
                phase = [one(summary, sample=sample, product=p, season=season, region=region) for p in products]
                text = "  ".join(f"{p} {r['phase_diff_h']:+.1f} h" for p, r in zip(products, phase) if r and np.isfinite(r["phase_diff_h"]))
                ax.set_title("amount phase vs gauge: " + text, loc="right", fontsize=5.4, color=INK_SECONDARY, pad=2)
            ax.set_xlim(0, 23)
            ax.set_xticks([0, 6, 12, 18, 23])
        label = ("All seasons" if stratum == "all" else stratum) if by == "season" else REGION_LABELS[stratum]
        axes[i, 0].text(-0.3, 0.5, label, transform=axes[i, 0].transAxes, rotation=90, ha="center", va="center", fontsize=6.8)
    for ax in axes[-1]:
        ax.set_xlabel("Hour of day (BJT, hour-ending)")
    handles = [Line2D([], [], color=INK, linewidth=1.8, label="Gauge")] + \
              [Line2D([], [], color=PRODUCT_COLORS[p], linewidth=1.5, label=p) for p in products]
    fig.legend(handles=handles, loc="upper left", ncol=len(handles), frameon=False, fontsize=6.8,
               bbox_to_anchor=(0.12, 1 - 0.62 / height))
    what = "season" if by == "season" else "landform region"
    fig.suptitle(f"Diurnal cycle by {what} — {SAMPLE_TITLES[sample]}", x=0.13, ha="left", fontsize=8.5, y=1 - 0.12 / height)
    note = "Hours with fewer than 100 cells are not drawn" + ("; FY4B has no 23:00 label." if sample == "all_products" else ".")
    fig.text(0.13, 1 - 0.32 / height, "Every source averaged over identical cells. Phase: first harmonic of the amount curve, "
             "+ = satellite later.\n" + note, ha="left", va="top", fontsize=6.0, color=INK_SECONDARY, linespacing=1.35)
    save(fig, f"fig6_diurnal_{by}_{sample}")


# --------------------------------------------------------------------------------------------------
# Figure 7: station correlation maps

def fig7_station_maps(stations: list[dict]) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 5.8))
    fig.subplots_adjust(left=0.12, right=0.97, top=0.88, bottom=0.06, wspace=0.08, hspace=0.18)
    cmap = ListedColormap(ORDINAL_BLUES + ["#0d366b"])
    bounds = [0.0, 0.2, 0.35, 0.5, 0.65, 1.0]
    norm = BoundaryNorm(bounds, cmap.N)
    landform = {s["station_id"]: s for s in read_rows(ROOT / "data" / "processed" / "covariates" / "station_landform.csv")}
    for i, sample in enumerate(SAMPLE_PRODUCTS):
        for j in range(3):
            ax = axes[i, j]
            products = SAMPLE_PRODUCTS[sample]
            if j >= len(products):
                ax.axis("off")
                continue
            product = products[j]
            rows = pick(stations, sample=sample, product=product, season="all")
            for region in REGIONS:
                members = [r for r in rows if r["region"] == region and np.isfinite(r["r_union_wet"])]
                lon = [landform[r["station_id"]]["lon"] for r in members]
                lat = [landform[r["station_id"]]["lat"] for r in members]
                ax.scatter(lon, lat, c=[r["r_union_wet"] for r in members], cmap=cmap, norm=norm, marker=REGION_MARKERS[region],
                           s=14, edgecolors=MUTED, linewidths=0.25, zorder=3)
            values = [r["r_union_wet"] for r in rows if np.isfinite(r["r_union_wet"])]
            ax.text(0.97, 0.03, f"median r {np.median(values):.2f}", transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=6, color=INK)
            ax.set_aspect(1 / np.cos(np.radians(32.2)), adjustable="datalim")
            ax.set_xticks([110, 111])
            ax.set_yticks([31.5, 32.5])
            ax.set_xticklabels(["110°E", "111°E"] if i == 1 or j >= 1 else [])
            ax.set_yticklabels(["31.5°N", "32.5°N"] if j == 0 else [])
            ax.tick_params(length=2, pad=1.5, labelsize=6)
            ax.set_title(product, loc="left", pad=3)
        axes[i, 0].text(-0.26, 0.5, SAMPLE_TITLES[sample].replace(" · ", "\n"), transform=axes[i, 0].transAxes, rotation=90,
                        ha="center", va="center", fontsize=6.5)
    # The empty sixth panel holds the colour scale and the marker key.
    box = axes[1, 2].get_position()
    cax = fig.add_axes([box.x0 + 0.02, box.y0 + 0.62 * box.height, box.width - 0.04, 0.014])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal", ticks=bounds)
    bar.outline.set_edgecolor(AXIS)
    bar.ax.tick_params(length=2, labelsize=6)
    bar.ax.set_title("Hourly Pearson r, hours where\ngauge or satellite ≥ 0.1 mm", fontsize=6.5, color=INK_SECONDARY,
                     loc="left", pad=4)
    handles = [Line2D([], [], linestyle="none", marker=REGION_MARKERS[r], markersize=4.5, markerfacecolor=MUTED,
                      markeredgecolor=MUTED, label=REGION_LABELS[r] + (" (11 gauges, not interpreted)" if r == "relief_ge1000" else ""))
               for r in REGIONS]
    fig.legend(handles=handles, loc="upper left", ncol=1, frameon=False, fontsize=6.3,
               bbox_to_anchor=(box.x0 + 0.01, box.y0 + 0.45 * box.height), title="Marker: landform region",
               title_fontsize=6.5, alignment="left")
    fig.suptitle("Hourly correlation at each gauge", x=0.08, ha="left", fontsize=9, y=0.975)
    fig.text(0.08, 0.93, "All seasons. Marker shape: landform region. Shared dry hours are excluded so they do not inflate r.",
             ha="left", fontsize=6.3, color=INK_SECONDARY)
    save(fig, "fig7_station_correlation_maps")


# --------------------------------------------------------------------------------------------------
# Figure 8: region x season scorecard

def fig8_scorecard(metrics, events, diurnal_summary, station_summary, sample: str) -> None:
    products = SAMPLE_PRODUCTS[sample]
    regions = ["all", "relief_lt500", "relief_500_1000"]

    def flagged(row, key):
        """(value, low_sample) from one summary row; a missing row is (None, True)."""
        return (None, True) if row is None else (row.get(key), bool(row.get("low_sample", False)))

    def event_row(p, s, g):
        return one(events, min_dry_gap_h=3.0, sample=sample, product=p, season=s, region=g, event_class="all")

    panels = [
        ("Relative bias, gauge-wet hours (%)", 0.0,
         lambda p, s, g: flagged(one(metrics, sample=sample, product=p, season=s, region=g, gauge_class="all_wet"), "RB_pct")),
        ("Heavy+ hours in the same class", 1.0,
         lambda p, s, g: (_pooled_hit(metrics, sample, p, s, g), False)),
        ("Event detected (POD)", 1.0, lambda p, s, g: flagged(event_row(p, s, g), "POD_event")),
        ("Event peak within ±1 h", 1.0, lambda p, s, g: flagged(event_row(p, s, g), "peak_within_1h")),
        ("Median station r (wet hours)", 1.0,
         lambda p, s, g: flagged(one(station_summary, sample=sample, product=p, season=s, region=g), "median_r_union_wet")),
        ("Diurnal phase error (h, + = late)", 0.0,
         lambda p, s, g: flagged(one(diurnal_summary, sample=sample, product=p, season=s, region=g), "phase_diff_h")),
    ]
    fig, axes = plt.subplots(len(panels), len(regions), figsize=(7.6, 9.6), sharey="row")
    fig.subplots_adjust(left=0.11, right=0.99, top=0.88, bottom=0.05, wspace=0.08, hspace=0.45)
    x = np.arange(len(SEASONS))
    offsets = np.linspace(-0.2, 0.2, len(products))
    for i, (title, ideal, getter) in enumerate(panels):
        for j, region in enumerate(regions):
            ax = axes[i, j]
            clean_axis(ax)
            ax.axhline(ideal, color=AXIS, linewidth=0.9, zorder=1)
            for offset, product in zip(offsets, products):
                pairs = [getter(product, season, region) for season in SEASONS]
                values = [np.nan if v is None else v for v, _ in pairs]
                dots(ax, x + offset, values, PRODUCT_COLORS[product], hollow=[low for _, low in pairs])
            ax.set_xticks(x)
            ax.set_xticklabels(["All", "MAM", "JJA", "SON", "DJF"] if i == len(panels) - 1 else [], fontsize=6)
            if i == 0:
                ax.set_title(REGION_LABELS[region], pad=4)
            if j == 0:
                ax.set_ylabel(title, fontsize=6.2)
    product_legend(fig, products, loc="upper right", bbox_to_anchor=(0.99, 0.985))
    fig.suptitle(f"Scorecard by season and landform region — {SAMPLE_TITLES[sample]}", x=0.11, ha="left", fontsize=9, y=0.985)
    fig.text(0.11, 0.962, "Relief ≥ 1000 m (11 gauges) is below the 20-gauge minimum and not shown. Region differences mix "
             "terrain with latitude:\nrugged gauges lie further south.", ha="left", va="top", fontsize=6.2,
             color=INK_SECONDARY, linespacing=1.35)
    fig.text(0.11, 0.934, CAVEAT, ha="left", va="top", fontsize=5.6, color=MUTED, linespacing=1.35)
    save(fig, f"fig8_scorecard_{sample}")


def _pooled_hit(metrics, sample, product, season, region):
    """Class-hit share pooled over the heavy, rainstorm and severe classes (weighted by n)."""
    rows = [one(metrics, sample=sample, product=product, season=season, region=region, gauge_class=c)
            for c in ("heavy", "rainstorm", "severe_rainstorm")]
    rows = [r for r in rows if r and r["n"] > 0]
    n = sum(r["n"] for r in rows)
    return sum(r["class_hit"] * r["n"] for r in rows) / n if n >= 100 else None


# --------------------------------------------------------------------------------------------------
# Figure 9: regional-mean series

SERIES_PANELS = [("r", "Correlation r", 1.0), ("KGE", "KGE", 1.0), ("RB_pct", "Relative bias (%)", 0.0)]


def fig9_regional_series(regional: list[dict]) -> None:
    regions = ["all", "relief_lt500", "relief_500_1000"]
    timescales = [1.0, 3.0, 24.0]
    fig, axes = plt.subplots(len(SERIES_PANELS), len(regions), figsize=(7.4, 5.6), sharey="row")
    fig.subplots_adjust(left=0.1, right=0.99, top=0.86, bottom=0.09, wspace=0.08, hspace=0.35)
    x = np.arange(len(timescales))
    for i, (key, title, ideal) in enumerate(SERIES_PANELS):
        for j, region in enumerate(regions):
            ax = axes[i, j]
            clean_axis(ax)
            ax.axhline(ideal, color=AXIS, linewidth=0.9, zorder=1)
            for sample, hollow in (("gpm_gsmap_full", False), ("all_products", True)):
                products = SAMPLE_PRODUCTS[sample]
                for product in products:
                    # Full-record dots on the left of each tick, FY4B-hours dots on the right, never sharing a position.
                    offset = {("gpm_gsmap_full", "GPM"): -0.3, ("gpm_gsmap_full", "GSMaP"): -0.17,
                              ("all_products", "FY4B"): 0.03, ("all_products", "GPM"): 0.16,
                              ("all_products", "GSMaP"): 0.29}[(sample, product)]
                    values = [(one(regional, sample=sample, product=product, region=region, season="all", timescale_h=t) or {}).get(key, np.nan)
                              for t in timescales]
                    dots(ax, x + offset, values, PRODUCT_COLORS[product], hollow=[hollow] * len(values))
            ax.set_xticks(x)
            ax.set_xlim(-0.5, len(timescales) - 0.5)
            ax.set_xticklabels(["1 h", "3 h", "24 h"] if i == len(SERIES_PANELS) - 1 else [])
            if i == 0:
                ax.set_title(REGION_LABELS[region], pad=4)
            if j == 0:
                ax.set_ylabel(title, fontsize=6.5)
    extra = [Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=MUTED, markeredgecolor=MUTED,
                    label="full record"),
             Line2D([], [], linestyle="none", marker="o", markersize=5, markerfacecolor=SURFACE, markeredgecolor=MUTED,
                    label="FY4B hours")]
    product_legend(fig, ["FY4B", "GPM", "GSMaP"], extra, loc="upper right", bbox_to_anchor=(0.99, 0.975))
    fig.suptitle("Regional-mean series across timescales", x=0.1, ha="left", fontsize=9, y=0.975)
    fig.text(0.1, 0.94, "Hourly means over each region's gauges (≥ 90% reporting), accumulated to met-day 3-h blocks and "
             "days, all seasons.\nFilled (left of each tick): full 2022–2024 record; hollow (right): the hours FY4B covers "
             "- no 24-h value, since FY4B never completes a day.", ha="left", va="top", fontsize=6.0,
             color=INK_SECONDARY, linespacing=1.35)
    save(fig, "fig9_regional_series")


def main() -> None:
    print(f"Drawing figures from {RUN_DIR}")
    metrics = read_rows(RUN_DIR / "intensity_class_metrics.csv")
    confusion = read_rows(RUN_DIR / "intensity_confusion.csv")
    false_alarms = read_rows(RUN_DIR / "false_alarm_metrics.csv")
    events = read_rows(RUN_DIR / "event_timing_summary.csv")
    diurnal = read_rows(RUN_DIR / "diurnal_cycle.csv")
    diurnal_scores = read_rows(RUN_DIR / "diurnal_summary.csv")
    stations = read_rows(RUN_DIR / "station_hourly_correlation.csv")
    station_summary = read_rows(RUN_DIR / "station_correlation_summary.csv")
    regional = read_rows(RUN_DIR / "regional_series_metrics.csv")

    fig1_landform()
    for sample in SAMPLE_PRODUCTS:
        fig2_intensity(metrics, sample)
    fig3_confusion(confusion)
    fig4_false_alarms(false_alarms)
    fig5_events(events)
    for sample in SAMPLE_PRODUCTS:
        fig6_diurnal(diurnal, diurnal_scores, sample, "season")
    fig6_diurnal(diurnal, diurnal_scores, "gpm_gsmap_full", "region")
    fig7_station_maps(stations)
    for sample in SAMPLE_PRODUCTS:
        fig8_scorecard(metrics, events, diurnal_scores, station_summary, sample)
    fig9_regional_series(regional)


if __name__ == "__main__":
    main()
