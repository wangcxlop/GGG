"""Draw the heavy-rainfall event evaluation figures.

Reads the CSVs written by `scripts/run_heavy_rain_event_evaluation.jl` and writes PNG (300 dpi)
and PDF figures to `output/heavy_rain_events/figures/`:

    py -3.13 scripts/plot_heavy_rain_events.py

The gauge-point maps (fig1) and station error maps (fig2) are the primary spatial comparison.
The IDW maps (fig3) are an auxiliary view of satellite values sampled at the gauge locations, not
native satellite fields; no metric is computed from them. Every panel of a kind shares one colour
scale, and every IDW surface was built with the same settings.

Colours follow a validated palette: a five-class single-hue blue ramp for precipitation, a
blue/red diverging ramp with a neutral grey centre for satellite-minus-gauge error (each arm
lightness-matched and checked as an ordinal ramp), and the first three categorical slots, which
validate all-pairs for colour-vision deficiency, for the three products.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "output" / "heavy_rain_events"
FIG_DIR = RUN_DIR / "figures"
PRODUCTS = ["FY4B", "GPM", "GSMaP"]
WINDOW = "common_hours"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
PRODUCT_COLORS = {"FY4B": "#2a78d6", "GPM": "#eb6834", "GSMaP": "#1baf7a"}
EVENT_MARKERS = ["o", "s", "^", "D", "v", "P"]

# CMA 24-h rain classes. No accumulation in these events reaches 250 mm, so 100-250 is the top class.
RAIN_BOUNDS = [0, 10, 25, 50, 100, 250]
RAIN_CMAP = ListedColormap(["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"])
RAIN_CMAP.set_over("#0d366b")
RAIN_NORM = BoundaryNorm(RAIN_BOUNDS, RAIN_CMAP.N)

# Satellite minus gauge: red = satellite drier, grey = within 10 mm, blue = satellite wetter.
ERROR_BOUNDS = [-250, -100, -50, -25, -10, 10, 25, 50, 100, 250]
ERROR_CMAP = ListedColormap([
    "#681114", "#9e3432", "#dd716a", "#f5938c", "#f0efec",
    "#86b6ef", "#5598e7", "#1c5cab", "#0d366b",
])
ERROR_CMAP.set_under("#681114")
ERROR_CMAP.set_over("#0d366b")
ERROR_NORM = BoundaryNorm(ERROR_BOUNDS, ERROR_CMAP.N)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK_SECONDARY,
    "axes.titlesize": 8,
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


def read_rows(name: str) -> list[dict[str, str]]:
    with open(RUN_DIR / name, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_events() -> list[dict]:
    """Selected events, widespread first, each group in date order."""
    events = [
        {"day": row["day"], "type": row["event_type"], "n_heavy": int(row["n_heavy"]),
         "max_mm": float(row["max_mm"]), "hours": int(row["common_hours"])}
        for row in read_rows("selected_events.csv")
    ]
    order = {"widespread": 0, "localized": 1}
    return sorted(events, key=lambda event: (order[event["type"]], event["day"]))


def load_stations() -> dict[tuple[str, str], dict[str, np.ndarray]]:
    """Common-window station totals keyed by (event day, product)."""
    columns: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in read_rows("station_event_totals.csv"):
        if row["window"] != WINDOW:
            continue
        entry = columns[(row["event_day"], row["product"])]
        for key in ("lon", "lat", "obs_mm", "sat_mm", "diff_mm"):
            entry[key].append(float(row[key]))
    return {key: {k: np.array(v) for k, v in value.items()} for key, value in columns.items()}


def load_metrics() -> dict[tuple[str, str, str], dict[str, float]]:
    metrics = {}
    for row in read_rows("event_metrics.csv"):
        values = {}
        for key, value in row.items():
            try:
                values[key] = float(value)
            except ValueError:
                values[key] = value
        metrics[(row["event_day"], row["product"], row["window"])] = values
    return metrics


def load_bounds() -> tuple[float, float, float, float]:
    row = read_rows("idw_settings.csv")[0]
    return tuple(float(row[key]) for key in ("west", "east", "south", "north"))


def event_label(event: dict) -> str:
    return (f"{event['day']} · {event['type']}\n"
            f"{event['n_heavy']} gauges > 50 mm · {event['hours']} h")


def save(fig: plt.Figure, name: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300)
    fig.savefig(FIG_DIR / f"{name}.pdf")
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


def style_map(ax: plt.Axes, bounds, show_x: bool, show_y: bool) -> None:
    west, east, south, north = bounds
    ax.set_xlim(west, east)
    ax.set_ylim(south, north)
    ax.set_aspect(1 / np.cos(np.radians((south + north) / 2)), adjustable="box")
    ax.set_xticks([110, 111])
    ax.set_yticks([32, 33])
    ax.set_xticklabels(["110°E", "111°E"] if show_x else [])
    ax.set_yticklabels(["32°N", "33°N"] if show_y else [])
    ax.tick_params(length=2, pad=1.5)


def signed(value: float) -> str:
    """Whole-number signed label that never prints "-0"."""
    return f"{round(value):+d}"


def annotate(ax: plt.Axes, text: str) -> None:
    """Metric box in the lower-right corner, which the gauge network leaves empty."""
    ax.text(0.97, 0.03, text, transform=ax.transAxes, ha="right", va="bottom", fontsize=5.6,
            color=INK, linespacing=1.15,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": SURFACE, "edgecolor": "none", "alpha": 0.85})


def row_label(ax: plt.Axes, event: dict) -> None:
    ax.text(-0.36, 0.5, event_label(event), transform=ax.transAxes, rotation=90,
            ha="center", va="center", fontsize=6.2, color=INK, linespacing=1.3)


def horizontal_colorbar(fig, rect, cmap, norm, bounds, label, extend="neither"):
    cax = fig.add_axes(rect)
    bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                       spacing="uniform", ticks=bounds, extend=extend)
    bar.outline.set_edgecolor(AXIS)
    bar.outline.set_linewidth(0.6)
    bar.ax.tick_params(length=2, labelsize=6.5)
    bar.set_label(label, fontsize=7, color=INK_SECONDARY)
    return bar


def points(ax, lon, lat, values, cmap, norm, size=9.0):
    """Gauge dots, drawn low-to-high magnitude so the heaviest totals stay on top."""
    order = np.argsort(np.abs(values))
    ax.scatter(lon[order], lat[order], c=values[order], cmap=cmap, norm=norm, s=size,
               edgecolors=MUTED, linewidths=0.2, zorder=3)


def fig1_point_maps(events, stations, metrics, bounds) -> None:
    columns = ["Gauge"] + PRODUCTS
    fig, axes = plt.subplots(len(events), 4, figsize=(7.0, 11.6))
    fig.subplots_adjust(left=0.13, right=0.985, top=0.935, bottom=0.075, wspace=0.06, hspace=0.08)
    for i, event in enumerate(events):
        gauge = stations[(event["day"], "GPM")]
        for j, source in enumerate(columns):
            ax = axes[i, j]
            style_map(ax, bounds, show_x=i == len(events) - 1, show_y=j == 0)
            if source == "Gauge":
                points(ax, gauge["lon"], gauge["lat"], gauge["obs_mm"], RAIN_CMAP, RAIN_NORM)
                annotate(ax, f"mean {gauge['obs_mm'].mean():.0f} mm\nmax {gauge['obs_mm'].max():.0f} mm")
            else:
                data = stations[(event["day"], source)]
                m = metrics[(event["day"], source, WINDOW)]
                points(ax, data["lon"], data["lat"], data["sat_mm"], RAIN_CMAP, RAIN_NORM)
                annotate(ax, f"r {m['r']:.2f}\nRMSE {m['RMSE']:.0f} mm\nBias {signed(m['Bias'])} mm")
            if i == 0:
                ax.set_title(source, pad=4)
        row_label(axes[i, 0], event)
    fig.suptitle("Heavy-rain events: gauge-observed vs satellite precipitation at the gauge locations",
                 fontsize=9, y=0.985)
    fig.text(0.5, 0.958, "Each dot is one gauge; satellite values are the pixels sampled at the same gauges. "
             "Totals over the hours all sources cover (08-08 BJT day).",
             ha="center", fontsize=6.5, color=INK_SECONDARY)
    horizontal_colorbar(fig, [0.30, 0.035, 0.52, 0.011], RAIN_CMAP, RAIN_NORM, RAIN_BOUNDS,
                        "Accumulated precipitation (mm), CMA daily rain classes")
    save(fig, "fig1_event_point_maps")


def fig2_error_maps(events, stations, metrics, bounds) -> None:
    fig, axes = plt.subplots(len(events), 3, figsize=(5.5, 11.6))
    fig.subplots_adjust(left=0.165, right=0.98, top=0.935, bottom=0.085, wspace=0.06, hspace=0.08)
    for i, event in enumerate(events):
        for j, product in enumerate(PRODUCTS):
            ax = axes[i, j]
            style_map(ax, bounds, show_x=i == len(events) - 1, show_y=j == 0)
            data = stations[(event["day"], product)]
            m = metrics[(event["day"], product, WINDOW)]
            points(ax, data["lon"], data["lat"], data["diff_mm"], ERROR_CMAP, ERROR_NORM)
            annotate(ax, f"Bias {signed(m['Bias'])} mm\nCRMSE {m['CRMSE']:.0f} mm")
            if i == 0:
                ax.set_title(f"{product} − gauge", pad=4)
        row_label(axes[i, 0], event)
    fig.suptitle("Station-level satellite-minus-gauge error", fontsize=9, y=0.985)
    fig.text(0.5, 0.958, "Same stations and hours as the point maps; one shared scale for every panel.",
             ha="center", fontsize=6.5, color=INK_SECONDARY)
    bar = horizontal_colorbar(fig, [0.20, 0.04, 0.72, 0.011], ERROR_CMAP, ERROR_NORM, ERROR_BOUNDS,
                              "Satellite − gauge (mm)", extend="both")
    bar.ax.text(0.0, 1.6, "← satellite drier", transform=bar.ax.transAxes, ha="left", fontsize=6.2, color=INK_SECONDARY)
    bar.ax.text(1.0, 1.6, "satellite wetter →", transform=bar.ax.transAxes, ha="right", fontsize=6.2, color=INK_SECONDARY)
    save(fig, "fig2_event_error_maps")


def load_surfaces() -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    columns: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in read_rows("idw_surfaces.csv"):
        entry = columns[(row["event_day"], row["source"])]
        entry["lon"].append(float(row["lon"]))
        entry["lat"].append(float(row["lat"]))
        entry["value"].append(float(row["value"]))
    surfaces = {}
    for key, entry in columns.items():
        lon, lat, value = (np.array(entry[k]) for k in ("lon", "lat", "value"))
        lons, lon_index = np.unique(np.round(lon, 6), return_inverse=True)
        lats, lat_index = np.unique(np.round(lat, 6), return_inverse=True)
        grid = np.full((lats.size, lons.size), np.nan)
        grid[lat_index, lon_index] = value
        surfaces[key] = (lons, lats, grid)
    return surfaces


# IDW extrapolates to the frame edge; beyond this distance from every gauge a surface would show
# colour where no sample exists, so those cells are left blank in every panel alike.
MAX_GAUGE_DISTANCE_KM = 20.0


def nearest_gauge_km(lons: np.ndarray, lats: np.ndarray, gauge_lon: np.ndarray, gauge_lat: np.ndarray) -> np.ndarray:
    """Great-circle distance from each grid cell `[lat, lon]` to its nearest gauge, in km."""
    lon_grid, lat_grid = np.meshgrid(np.radians(lons), np.radians(lats))
    glon, glat = np.radians(gauge_lon), np.radians(gauge_lat)
    a = (np.sin((lat_grid[..., None] - glat) / 2) ** 2
         + np.cos(lat_grid[..., None]) * np.cos(glat) * np.sin((lon_grid[..., None] - glon) / 2) ** 2)
    return (2 * 6378.388 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))).min(axis=-1)


def fig3_idw_maps(events, stations, bounds) -> None:
    settings = read_rows("idw_settings.csv")[0]
    surfaces = load_surfaces()
    columns = ["Gauge"] + PRODUCTS
    fig, axes = plt.subplots(len(events), 4, figsize=(7.0, 11.6))
    fig.subplots_adjust(left=0.13, right=0.985, top=0.905, bottom=0.075, wspace=0.06, hspace=0.08)
    for i, event in enumerate(events):
        gauge = stations[(event["day"], "GPM")]
        lons, lats, _ = surfaces[(event["day"], "Gauge")]
        far = nearest_gauge_km(lons, lats, gauge["lon"], gauge["lat"]) > MAX_GAUGE_DISTANCE_KM
        for j, source in enumerate(columns):
            ax = axes[i, j]
            style_map(ax, bounds, show_x=i == len(events) - 1, show_y=j == 0)
            lons, lats, grid = surfaces[(event["day"], source)]
            ax.pcolormesh(lons, lats, np.ma.masked_where(far | np.isnan(grid), grid), cmap=RAIN_CMAP,
                          norm=RAIN_NORM, shading="nearest", rasterized=True, zorder=1)
            values = gauge["obs_mm"] if source == "Gauge" else stations[(event["day"], source)]["sat_mm"]
            order = np.argsort(values)
            ax.scatter(gauge["lon"][order], gauge["lat"][order], c=values[order], cmap=RAIN_CMAP,
                       norm=RAIN_NORM, s=3.5, edgecolors=INK, linewidths=0.15, zorder=3)
            if i == 0:
                ax.set_title(source, pad=4)
        row_label(axes[i, 0], event)
    fig.suptitle("Auxiliary view: IDW visualization of satellite values sampled at gauge locations",
                 fontsize=9, y=0.985)
    fig.text(0.5, 0.961, "Not native satellite fields: every surface is interpolated from the same gauge locations, "
             "so it cannot show structure between gauges.\nMetrics are not computed from these surfaces. "
             f"Identical settings for every panel: inverse distance weighting, power {float(settings['power']):g}, "
             f"{settings['neighbors']},\n{float(settings['step_deg']):g}° grid, cells over "
             f"{MAX_GAUGE_DISTANCE_KM:g} km from any gauge left blank. Dots are the sampled values.",
             ha="center", va="top", fontsize=6.3, color=INK_SECONDARY, linespacing=1.35)
    horizontal_colorbar(fig, [0.30, 0.035, 0.52, 0.011], RAIN_CMAP, RAIN_NORM, RAIN_BOUNDS,
                        "Accumulated precipitation (mm), CMA daily rain classes")
    save(fig, "fig3_event_idw_maps")


def fig4_scatter(events, stations, metrics) -> None:
    fig, axes = plt.subplots(len(events), 3, figsize=(5.6, 11.6))
    fig.subplots_adjust(left=0.2, right=0.98, top=0.925, bottom=0.05, wspace=0.18, hspace=0.34)
    for i, event in enumerate(events):
        row = [stations[(event["day"], product)] for product in PRODUCTS]
        limit = 1.05 * max(max(d["obs_mm"].max(), d["sat_mm"].max()) for d in row)
        for j, product in enumerate(PRODUCTS):
            ax = axes[i, j]
            data = row[j]
            m = metrics[(event["day"], product, WINDOW)]
            ax.grid(color=GRID, linewidth=0.5, zorder=0)
            ax.plot([0, limit], [0, limit], color=AXIS, linewidth=0.8, zorder=1)
            ax.scatter(data["obs_mm"], data["sat_mm"], s=8, color=PRODUCT_COLORS[product], alpha=0.8,
                       edgecolors=SURFACE, linewidths=0.3, zorder=2)
            ax.set_xlim(0, limit)
            ax.set_ylim(0, limit)
            ax.set_aspect("equal", adjustable="box")
            ax.tick_params(length=2, pad=1.5, labelsize=6)
            for spine in ("top", "right"):
                ax.spines[spine].set_visible(False)
            # Above the panel rather than inside it: no corner of a scatter is reliably empty.
            ax.set_title(f"r {m['r']:.2f}  ·  RMSE {m['RMSE']:.0f} mm\nBias {signed(m['Bias'])} mm  ·  n {int(m['n'])}",
                         loc="left", fontsize=5.8, color=INK_SECONDARY, pad=3, linespacing=1.2)
            if i == 0:
                ax.text(0.5, 1.3, product, transform=ax.transAxes, ha="center", va="bottom", fontsize=8)
            if i == len(events) - 1:
                ax.set_xlabel("Gauge (mm)")
            if j == 0:
                ax.set_ylabel("Satellite (mm)")
        axes[i, 0].text(-0.62, 0.5, event_label(event), transform=axes[i, 0].transAxes, rotation=90,
                        ha="center", va="center", fontsize=6.2, color=INK, linespacing=1.3)
    fig.suptitle("Satellite vs gauge accumulations at the same gauges (1:1 line)", fontsize=9, y=0.985)
    save(fig, "fig4_event_scatter")


METRIC_PANELS = [
    ("r", "Spatial correlation r", 1.0),
    ("rho", "Spearman ρ", 1.0),
    ("KGE", "KGE", 1.0),
    ("CSI_50", "CSI at 50 mm", 1.0),
    ("RMSE", "RMSE (mm)", 0.0),
    ("RB_pct", "Relative bias (%)", 0.0),
    ("cv_ratio", "Spatial CV ratio (satellite / gauge)", 1.0),
    ("centroid_shift_km", "Rain-centre shift (km)", 0.0),
]


def fig5_metric_summary(events, metrics) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(10.0, 5.4))
    fig.subplots_adjust(left=0.055, right=0.99, top=0.84, bottom=0.14, wspace=0.28, hspace=0.42)
    x = np.arange(len(events))
    offsets = {"FY4B": -0.2, "GPM": 0.0, "GSMaP": 0.2}
    n_widespread = sum(event["type"] == "widespread" for event in events)
    tick_labels = [f"{event['day'][5:]}\n{event['day'][:4]}" for event in events]
    for ax, (key, title, ideal) in zip(axes.flat, METRIC_PANELS):
        ax.grid(axis="y", color=GRID, linewidth=0.5, zorder=0)
        ax.axhline(ideal, color=AXIS, linewidth=0.9, zorder=1)
        ax.axvline(n_widespread - 0.5, color=GRID, linewidth=0.8, zorder=0)
        for product in PRODUCTS:
            values = [metrics[(event["day"], product, WINDOW)][key] for event in events]
            ax.scatter(x + offsets[product], values, s=26, color=PRODUCT_COLORS[product],
                       edgecolors=SURFACE, linewidths=0.8, zorder=3)
        ax.set_title(title, loc="left", fontsize=7.5, pad=13)
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, fontsize=6)
        ax.set_xlim(-0.6, len(events) - 0.4)
        ax.tick_params(length=2, pad=1.5, labelsize=6)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.text((n_widespread - 1) / 2, 1.0, "widespread", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=5.8, color=MUTED)
        ax.text(n_widespread + (len(events) - n_widespread - 1) / 2, 1.0, "localized",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=5.8, color=MUTED)
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=5.5, markerfacecolor=PRODUCT_COLORS[p],
                      markeredgecolor=SURFACE, label=p) for p in PRODUCTS]
    fig.legend(handles=handles, loc="upper right", ncol=3, frameon=False, fontsize=7, bbox_to_anchor=(0.99, 0.975))
    fig.suptitle("Event-scale spatial metrics by product", x=0.055, ha="left", fontsize=9, y=0.965)
    fig.text(0.055, 0.915, "Gauge-pixel pairs over each event's common hours. The solid grey line marks "
             "the perfect score.", ha="left", fontsize=6.5, color=INK_SECONDARY)
    save(fig, "fig5_metric_summary")


def fig6_taylor(events, metrics) -> None:
    r_min = -0.3
    theta_max = np.degrees(np.arccos(r_min))
    sd_values = [metrics[(e["day"], p, WINDOW)]["sd_ratio"] for e in events for p in PRODUCTS]
    r_max = max(1.5, np.ceil(max(sd_values) * 2) / 2 + 0.25)

    fig = plt.figure(figsize=(7.0, 5.2))
    ax = fig.add_axes([0.06, 0.1, 0.62, 0.78], projection="polar")
    ax.set_thetamin(0)
    ax.set_thetamax(theta_max)
    ax.set_rlim(0, r_max)
    ax.grid(color=GRID, linewidth=0.5)
    ax.spines["polar"].set_edgecolor(AXIS)
    r_ticks = [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99]
    ax.set_thetagrids(np.degrees(np.arccos(r_ticks)), [f"{r:g}" for r in r_ticks], fontsize=6.5, color=INK_SECONDARY)
    ax.set_rgrids(np.arange(0.5, r_max, 0.5), fontsize=6.5, color=INK_SECONDARY, angle=theta_max)
    ax.tick_params(pad=2)

    # Centred-RMSE contours around the gauge reference point (r = 1, normalized std = 1).
    phi = np.linspace(0, np.pi, 400)
    for level in (0.5, 1.0, 1.5):
        xs, ys = 1 + level * np.cos(phi), level * np.sin(phi)
        radius, angle = np.hypot(xs, ys), np.arctan2(ys, xs)
        keep = (radius <= r_max) & (angle <= np.radians(theta_max))
        ax.plot(angle[keep], radius[keep], color=AXIS, linewidth=0.6, zorder=1)
        # Labelled near r = 1, where no event lands, so no label sits under a marker.
        label_at = np.argmin(np.abs(np.degrees(angle[keep]) - 4)) if keep.any() else None
        if label_at is not None:
            ax.text(angle[keep][label_at], radius[keep][label_at], f"CRMSE {level:g}", fontsize=5.5,
                    color=MUTED, ha="center", va="center",
                    bbox={"boxstyle": "round,pad=0.1", "facecolor": SURFACE, "edgecolor": "none"})

    ax.scatter([0], [1], marker="*", s=90, color=INK, zorder=5, clip_on=False)
    ax.text(0, 1.0, "  Gauge", fontsize=6.5, color=INK, ha="left", va="bottom")
    for event, marker in zip(events, EVENT_MARKERS):
        for product in PRODUCTS:
            m = metrics[(event["day"], product, WINDOW)]
            ax.scatter(np.arccos(np.clip(m["r"], -1, 1)), m["sd_ratio"], marker=marker, s=34,
                       color=PRODUCT_COLORS[product], edgecolors=SURFACE, linewidths=0.8, zorder=3)
    ax.text(np.radians(theta_max / 2), r_max * 1.13, "Spatial correlation", rotation=theta_max / 2 - 90,
            ha="center", va="center", fontsize=7, color=INK_SECONDARY)
    ax.text(0, r_max / 2, "\n\nNormalized standard deviation (satellite / gauge)", ha="center", va="top",
            fontsize=7, color=INK_SECONDARY)

    product_handles = [Line2D([], [], linestyle="none", marker="o", markersize=6, markerfacecolor=PRODUCT_COLORS[p],
                              markeredgecolor=SURFACE, label=p) for p in PRODUCTS]
    event_handles = [Line2D([], [], linestyle="none", marker=marker, markersize=5.5, markerfacecolor="none",
                            markeredgecolor=INK, label=f"{e['day']} ({e['type']})")
                     for e, marker in zip(events, EVENT_MARKERS)]
    fig.legend(handles=product_handles, title="Product", loc="upper left", bbox_to_anchor=(0.72, 0.8),
               frameon=False, fontsize=7, title_fontsize=7.5, alignment="left")
    fig.legend(handles=event_handles, title="Event", loc="upper left", bbox_to_anchor=(0.72, 0.56),
               frameon=False, fontsize=7, title_fontsize=7.5, alignment="left")
    fig.suptitle("Taylor diagram of event-scale spatial patterns", x=0.06, ha="left", fontsize=9, y=0.965)
    fig.text(0.06, 0.915, "Distance from the gauge star is the centred RMSE (normalized): pattern error with "
             "the mean bias removed.", ha="left", fontsize=6.5, color=INK_SECONDARY)
    save(fig, "fig6_taylor_diagram")


def main() -> None:
    events = load_events()
    stations = load_stations()
    metrics = load_metrics()
    bounds = load_bounds()
    print(f"Drawing {len(events)} events from {RUN_DIR}")
    fig1_point_maps(events, stations, metrics, bounds)
    fig2_error_maps(events, stations, metrics, bounds)
    fig3_idw_maps(events, stations, bounds)
    fig4_scatter(events, stations, metrics)
    fig5_metric_summary(events, metrics)
    fig6_taylor(events, metrics)


if __name__ == "__main__":
    main()
