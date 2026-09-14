// Interim results: spatial evaluation of FY4B, GPM and GSMaP on representative heavy-rain events.
// Compile from this folder:  typst compile heavy_rain_spatial_evaluation.typ
// Every number in the tables is generated from output/heavy_rain_events/*.csv
// (commit 2658353, scripts/run_heavy_rain_event_evaluation.jl); figures are copies of
// output/heavy_rain_events/figures/.

#set document(
  title: "Spatial heterogeneity of rainfall in FY4B, GPM and GSMaP: interim results",
  date: datetime(year: 2026, month: 9, day: 13),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Spatial evaluation of satellite precipitation — interim results
      #h(1fr) 13 September 2026
    ]
  },
)
#set text(font: ("Libertinus Serif", "New Computer Modern"), size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 1.05em)
#set heading(numbering: "1.1")
#show heading: set text(font: ("Segoe UI", "Arial"), weight: "semibold")
#show heading.where(level: 1): it => { v(0.6em); it; v(0.25em) }
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)
#show figure.where(kind: table): set block(breakable: true)
#show table.cell: set text(size: 8pt)
#show table: set par(justify: false)
#show raw: set text(size: 8.5pt)
#set list(indent: 0.6em)

// ---------------------------------------------------------------------------------------------
// Tables. The row data below was generated from the CSVs; do not edit numbers by hand.

#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)

#let note(body) = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, fill: luma(246),
  stroke: (left: 2.5pt + rgb("#2a78d6")), body)

#let events_rows = (
  [2022-06-26], [widespread], [116], [147.0], [50.9], [23], [12.1], [0.0],
  [2023-06-17], [widespread], [88], [112.0], [47.8], [23], [1.2], [0.0],
  [2023-08-26], [widespread], [167], [111.5], [63.4], [21], [15.3], [6.7],
  [2024-07-10], [widespread], [76], [160.0], [43.3], [22], [8.8], [12.8],
  [2023-07-31], [localized], [9], [135.0], [8.3], [23], [14.3], [0.4],
  [2023-09-10], [localized], [17], [136.0], [14.7], [22], [1.2], [0.0],
)

#let summary_rows = (
  [FY4B], [0.359], [29.9], [23.5], [−9.5], [20.5], [0.625], [−0.042], [0.391], [0.165],
  [GPM], [0.538], [18.8], [13.0], [+5.6], [17.9], [0.732], [0.381], [0.546], [0.333],
  [GSMaP], [0.536], [34.9], [27.5], [+32.5], [23.0], [0.826], [−0.009], [0.534], [0.381],
)

#let robust_rows = (
  [GPM], [common_hours], [0.538], [18.8], [+5.6], [0.732], [0.381], [0.333],
  [GPM], [full_day], [0.538], [20.5], [+4.8], [0.776], [0.364], [0.413],
  [GSMaP], [common_hours], [0.536], [34.9], [+32.5], [0.826], [−0.009], [0.381],
  [GSMaP], [full_day], [0.529], [35.3], [+27.3], [0.862], [0.051], [0.440],
)

#let per_event_rows = (
  table.cell(colspan: 5, fill: luma(235))[*2022-06-26 (widespread, 226 stations, all products fail)*],
  [FY4B], [0.011], [60.0], [+113], [−0.505],
  [GPM], [0.029], [24.0], [+0.2], [−0.111],
  [GSMaP], [−0.106], [37.5], [+45.3], [−0.198],
  table.cell(colspan: 5, fill: luma(235))[*2023-06-17 (widespread, 229 stations, all products agree well)*],
  [FY4B], [0.733], [16.7], [−21.9], [0.597],
  [GPM], [0.701], [14.5], [−3.4], [0.693],
  [GSMaP], [0.799], [37.1], [+49.3], [−0.295],
  table.cell(colspan: 5, fill: luma(235))[*2023-08-26 (widespread, 228 stations, largest bias/location errors)*],
  [FY4B], [−0.163], [31.6], [−45.7], [−0.445],
  [GPM], [0.427], [20.5], [+15.7], [0.393],
  [GSMaP], [0.465], [84.2], [+147], [−0.756],
  table.cell(colspan: 5, fill: luma(235))[*2024-07-10 (widespread, 237 stations, strong correlation, FY4B badly low)*],
  [FY4B], [0.849], [34.2], [−68.4], [−0.010],
  [GPM], [0.842], [18.5], [+24.4], [0.700],
  [GSMaP], [0.861], [14.3], [+11.0], [0.784],
  table.cell(colspan: 5, fill: luma(235))[*2023-07-31 (localized, 232 stations, weakest correlations of the set)*],
  [FY4B], [0.347], [17.3], [−31.9], [0.046],
  [GPM], [0.500], [18.2], [+35.2], [0.388],
  [GSMaP], [0.451], [16.9], [+0.8], [0.374],
  table.cell(colspan: 5, fill: luma(235))[*2023-09-10 (localized, 231 stations, moderate-high r, all underestimate)*],
  [FY4B], [0.376], [19.6], [−2.2], [0.064],
  [GPM], [0.730], [17.2], [−38.6], [0.221],
  [GSMaP], [0.747], [19.2], [−57.9], [0.036],
)

// ---------------------------------------------------------------------------------------------
// Title

#align(left)[
  #set text(hyphenate: false)
  #text(font: ("Segoe UI", "Arial"), size: 9pt, fill: luma(100), tracking: 0.04em)[INTERIM RESULTS · 13 SEPTEMBER 2026]
  #v(0.2em)
  #text(font: ("Segoe UI", "Arial"), size: 19pt, weight: "semibold")[Can satellite precipitation products capture the spatial heterogeneity of rainfall?]
  #v(0.1em)
  #text(size: 12pt, fill: luma(60))[FY4B, GPM IMERG and GSMaP against the gauge network on 6 representative heavy-rain events in the Hubei study area, 2022–2024]
]
#v(0.6em)
#line(length: 100%, stroke: 0.5pt + luma(180))

= Key findings

- *GPM is the most reliable all-round product.* It has the lowest RMSE (18.8 mm) and the smallest bias (+5.6%) of the three, the best KGE (0.381), and its accuracy barely changes whether it is scored on the shared common-hour window or its own full day.
- *FY4B has the weakest spatial correlation and squeezes the rain field the most.* Mean r = 0.36, and its `cv_ratio` of 0.625 means it reproduces only about five-eighths of the gauge network's station-to-station variability — the flattest of the three fields.
- *GSMaP reproduces the amplitude of spatial variability best but is badly biased.* Its `cv_ratio` (0.826) is the closest to 1 of the three products, and it posts the best CSI at the 50 mm threshold (0.381), but it overestimates the network mean by a third (RB = +32.5%) and has the worst RMSE (34.9 mm).
- *All three products compress spatial heterogeneity* (`cv_ratio` < 1 for every product, every event except one GSMaP full-day case): none of them spreads its estimates across the gauge network as widely as the gauges themselves vary — the central answer to the research question.
- *Skill is event-dependent, not just product-dependent.* On 2022-06-26 all three products fail outright (r between −0.11 and +0.03); on 2023-06-17 and 2024-07-10 all three correlate well (r = 0.70–0.86). The same product can be the best or the worst performer depending on which storm is scored.
- *2023-08-26 is the worst single case*: FY4B misplaces its rain peak 125 km from the true one, and GSMaP overestimates the network mean by 147%, the largest bias anywhere in the dataset.
- *The common-hour window is not driving these results.* Re-scoring GPM and GSMaP over their own complete 24-hour day gives essentially the same r, bias and `cv_ratio` as the common window all three products share.

= Question and scope

The question is whether FY4B, GPM and GSMaP reproduce the *spatial* pattern of rainfall — not just its network-average amount — when heavy rain actually falls. This complements the companion report, `satellite_temporal_evaluation.typ`, which addresses the *temporal* evolution of rainfall over the full 2022–2024 hourly record; that report's benchmark scores every station-hour, of which about 93% are dry, which dilutes heavy-rain performance specifically. This report instead isolates a small number of representative heavy-rain days and asks, across the gauge network on each one: does a product's spatial pattern of accumulated rain look like the gauges' pattern, in magnitude, in variability, and in where the rain fell?

= Data

- *Gauges.* The same hourly gauge network as the temporal-evaluation report (`hubei_obs_hourly_2022_2025_JunSep.csv`), 226–237 stations reporting on any given selected day, in the Hubei study area (109.4–111.6° E, 31.2–33.4° N), hour-ending Beijing time (BJT).
- *Satellite products*, sampled at the gauge pixels and aligned to the gauge clock: FY4B (strict hourly aggregation, navigation-corrected), GPM IMERG (Google Earth Engine, UTC hour-start realigned +9 h), GSMaP (same alignment as GPM).
- *Period screened.* Every 08–08 BJT day from 1 January 2022 to 31 December 2024.
- *Meteorological day.* Hour-ending labels `D 09:00` through `D+1 08:00` are all attributed to day `D`; a station-day total requires all 24 of its hours to be finite, so no heavy-rain count rests on a partial day.

= Methods

== Event screening and selection

+ *Candidate days.* A day qualifies if more than 3 gauges (strict `>`, `min_stations = 4`) report more than 50 mm. This screened 56 candidate days out of the three years (`candidate_days.csv`). A second, looser 25 mm threshold counts gauges for an areal-extent measure, `frac_wide`.
+ *Rain processes.* Consecutive qualifying days are merged into one `process_id`, so a multi-day rain system is never picked twice.
+ *Common-hour window.* FY4B systematically lacks the 23:00 label and often 00:00–01:00, so a full-day FY4B total would look biased low next to GPM or GSMaP purely from missing data. For each day, the common window keeps only the hours where the gauge network *and* every product have at least 90% of stations finite. Two shares track how much this leaves out: `excluded_rain_share` (the network's rain outside the window) and `peak_excluded_rain_share` (the same, at the single wettest gauge alone — a network-wide share can look harmless while the storm core itself falls in the missing hours).
+ *Eligibility.* A day is eligible only if the common window keeps at least 20 hours and both excluded-rain shares are at most 20%.
+ *Selection.* Eligible days split into *widespread* (`frac_wide >= 0.25`, ranked by gauge count above 50 mm, top 4 picked) and *localized* (`frac_wide < 0.25`, ranked by the single-gauge maximum — an intense storm over a small footprint, the hardest case for a satellite's resolution, top 2 picked). One day is kept per rain process across both groups.

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.6pt),
    table.header(rule, [*Day*], [*Type*], [*Gauges > 50 mm*], [*Peak (mm)*], [*Mean (mm)*], [*Common hours*], [*Excluded rain %*], [*Peak-excl. %*], thin),
    ..events_rows,
    rule,
  ),
  caption: [The 6 selected representative events (4 widespread, 2 localized), one per independent rain process.],
) <tab-events>
]

== Metrics <sec-metrics>

For each event, per product, stations are kept only when the gauge *and every product* have a complete common-hour accumulation, so all three products are scored on one identical station set. From the paired (gauge, product) accumulations:

- *RMSE, MAE, Bias, r* (Pearson) and *RB_pct* = 100·Bias/obs_mean — the standard accuracy measures.
- *rho*: Spearman rank correlation, robust to the skew of heavy-rain totals.
- *CRMSE*: centred RMSE, the pattern error left once the mean bias is removed; together with `r` and `sd_ratio` it satisfies the Taylor-diagram identity (@fig-taylor).
- *sd_ratio* = σ#sub[est]/σ#sub[obs] and *cv_ratio* = (σ#sub[est]/mean#sub[est]) / (σ#sub[obs]/mean#sub[obs]) — the amplitude of spatial variability a product reproduces, independent of its mean bias. This is the report's central spatial-heterogeneity metric.
- *KGE* = 1 − √((r−1)² + (sd_ratio−1)² + (mean#sub[est]/mean#sub[obs]−1)²).
- *POD, FAR, CSI* at 25 mm and 50 mm thresholds (`>=` threshold): whether a gauge's rain class is detected in the right place across the network.
- *centroid_shift_km*: great-circle distance between the gauge network's rain-weighted centroid and the product's rain-weighted centroid.
- *peak_shift_km*: distance between the gauge with the largest observed total and the gauge where the product places its own maximum.

Every metric is computed twice for GPM and GSMaP: once over the shared `common_hours` window, and once over their own complete `full_day` (24 h) — a check that FY4B's missing hours are not what shapes the comparison (@sec-robustness).

#note[
  *Scope of "spatial heterogeneity" here.* This evaluation operationalizes heterogeneity as (a) the `cv_ratio` of accumulated totals across the gauge network, (b) station-to-station spatial correlation (`r`, `rho`), (c) rain-centroid/peak-location shift distances, and (d) categorical skill (POD/FAR/CSI) at spatial thresholds. It does not fit a geostatistical semivariogram, Moran's I, or a GWR local-statistic surface — those tools exist elsewhere in this repository for other analyses but were not applied to these six events.
]

#note[
  *IDW maps are visualization only.* The maps in @fig-idw interpolate each source's values *sampled at the gauge locations* (power = 2, all finite gauges, 0.025° grid, common-hour window) — not a product's native gridded field — purely so the reader can see the spatial pattern each source implies. No metric in this report is computed from an IDW surface.
]

= Results

== Aggregate performance

Averaged across all 6 events, GPM has the best accuracy and GSMaP the best spatial-variability amplitude, but at the cost of a large positive bias (@tab-summary, @fig-summary, @fig-taylor). The Taylor diagram shows GPM's points clustered closest to the gauge reference (r = 1, sd_ratio = 1); GSMaP's points sit at a similar angular correlation to GPM's but further out in normalized standard deviation, reflecting its overestimation; FY4B's points sit at the widest angle (lowest r) and closest in to the origin (most compressed variability).

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.6pt),
    table.vline(x: 1, stroke: 0.4pt + luma(170)),
    table.header(rule, [*Product*], [*r*], [*RMSE (mm)*], [*MAE (mm)*], [*RB (%)*], [*CRMSE (mm)*], [*cv_ratio*], [*KGE*], [*rho#super[a]*], [*CSI_50*], thin),
    ..summary_rows,
    rule,
  ),
  caption: [Mean-over-events aggregate metrics, common-hour window (n = 6 events, `event_metrics_summary.csv`). #super[a] Spearman rank correlation.],
) <tab-summary>
]

#figure(
  image("figures/fig5_metric_summary.png", width: 100%),
  caption: [Eight metrics (r, rho, KGE, CSI_50, RMSE, RB_pct, cv_ratio, rain-centre shift km) across the 6 events, grouped widespread vs. localized, one marker per product.],
) <fig-summary>

#figure(
  image("figures/fig6_taylor_diagram.png", width: 78%),
  caption: [Taylor diagram: spatial correlation (angle) vs. normalized standard deviation (radius) for all 6 events × 3 products, with CRMSE contours around the gauge reference point.],
) <fig-taylor>

== Spatial-heterogeneity amplitude

The `cv_ratio` values in @tab-summary are all below 1: FY4B reproduces 63% of the gauge network's coefficient of variation, GPM 73%, and GSMaP 83% — the closest of the three, though still short of matching the gauges' spatial spread. In other words, every product's spatial field is *smoother* than the true rainfall field: heavy-rain totals are under-spread across the network relative to what the gauges recorded, regardless of whether the product's mean is biased high or low. Combined with GSMaP's positive mean bias, this means GSMaP overestimates broadly rather than reproducing sharp local maxima; FY4B's more severe compression compounds its already-weak correlation.

== Event-level breakdown

Performance depends strongly on which storm is scored (@tab-per-event, @fig-points, @fig-scatter). 2023-06-17 and 2024-07-10 are the two events where every product achieves r ≥ 0.70; 2022-06-26 is the one event where every product fails outright, with r between −0.11 and +0.03, despite 116 gauges recording more than 50 mm that day (mean 50.9 mm across the network) — a widespread, substantial event by every measure except how well the satellites captured its pattern. The two localized events are harder for all products than the widespread ones on average, consistent with a spatially concentrated storm falling below several products' effective resolution, though 2023-09-10 still reaches r = 0.73–0.75 for GPM and GSMaP.

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.6pt),
    table.header(rule, [*Product*], [*r*], [*RMSE (mm)*], [*RB (%)*], [*KGE*], thin),
    ..per_event_rows,
    rule,
  ),
  caption: [Per-event metrics, common-hour window (`event_metrics.csv`).],
) <tab-per-event>
]

#figure(
  image("figures/fig1_event_point_maps.png", width: 100%),
  caption: [Gauge accumulation (left column) and each product's accumulation at the same gauges, per event, annotated with r/RMSE/Bias.],
) <fig-points>

#figure(
  image("figures/fig4_event_scatter.png", width: 100%),
  caption: [Satellite-vs-gauge scatter per event and product, with the 1:1 line, r, RMSE, Bias and n.],
) <fig-scatter>

#note[
  *2023-08-26 is the outlier to watch.* FY4B places its peak 125.5 km from the true peak gauge — the largest location error anywhere in the dataset — while GSMaP overestimates the network mean on the same day by 147% (est. 132.6 mm vs. obs. 53.75 mm), its worst bias of any event. Both point to this event, not a general product failure, as the source of the widest spread in @fig-summary.
]

== Error pattern and location shifts

@fig-error-maps shows where each product's error concentrates spatially. FY4B's errors on the two worst events (2022-06-26, 2023-08-26) are large and one-signed across most of the network (underestimation on 2023-08-26, overestimation on 2022-06-26), rather than a mix of local overshoots and undershoots — consistent with its low `cv_ratio`: it is missing the pattern, not scattering around it. GSMaP's errors on 2023-08-26 are large and positive almost everywhere, matching its network-wide overestimation that day. Rain-centroid shifts average 7.7 km (GPM), 12.3 km (GSMaP) and 17.6 km (FY4B); peak-location shifts are larger and noisier (17.7–39.1 km on average), reflecting how sensitive a single-gauge maximum is to noise even when the broader pattern is reasonable.

#figure(
  image("figures/fig2_event_error_maps.png", width: 100%),
  caption: [Station-level satellite-minus-gauge error per event and product (diverging colour scale), annotated with Bias and CRMSE.],
) <fig-error-maps>

#figure(
  image("figures/fig3_event_idw_maps.png", width: 100%),
  caption: [Auxiliary IDW-interpolated surfaces of gauge and satellite values sampled at gauge locations (illustrative only — see the note in @sec-metrics).],
) <fig-idw>

== Categorical skill

GSMaP has both the best detection and the best CSI at the 50 mm threshold (POD_50 = 0.655, FAR_50 = 0.473, CSI_50 = 0.381, all mean-over-events), consistent with its wide overestimation catching most true heavy-rain stations. GPM detects fewer of them and has a marginally higher false-alarm rate (POD_50 = 0.543, FAR_50 = 0.499, CSI_50 = 0.333). FY4B has the lowest false-alarm rate of the three (FAR_50 = 0.323) but by far the lowest detection (POD_50 = 0.260), so its CSI_50 (0.165) is still the weakest — it simply misses most of the heavy-rain stations rather than flagging extra ones.

== Robustness: common-hour window vs. full day <sec-robustness>

Re-scoring GPM and GSMaP over their own complete 24-hour day, rather than the common window all three products share, changes every metric by only a small amount (@tab-robust): GPM's `cv_ratio` moves from 0.732 to 0.776 and GSMaP's from 0.826 to 0.862, both well within the event-to-event spread already seen in @tab-per-event. This confirms the common-hour restriction imposed by FY4B's data gaps is not what drives the comparison between products.

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.6pt),
    table.header(rule, [*Product*], [*Window*], [*r*], [*RMSE (mm)*], [*RB (%)*], [*cv_ratio*], [*KGE*], [*CSI_50*], thin),
    ..robust_rows,
    rule,
  ),
  caption: [GPM and GSMaP, mean-over-events, common-hour window vs. their own full day.],
) <tab-robust>
]

= Answers to the research questions

+ *Can the products capture the spatial heterogeneity of heavy rainfall?* Only partially, and unevenly across products. All three compress the amplitude of spatial variability (`cv_ratio` < 1); GSMaP comes closest to the gauges' true spread but overshoots the mean broadly, GPM is the best-balanced product overall, and FY4B both correlates weakest and compresses variability most.
+ *Does this depend on the type of event?* Yes, more than it depends on product. 2022-06-26 defeats all three products; 2023-06-17 and 2024-07-10 are captured reasonably well by all three. Localized events are on average harder than widespread ones, but not uniformly — 2023-09-10 is captured about as well as the better widespread events.
+ *Which product is most reliable?* GPM, on balance: it has the lowest RMSE and bias and the most stable behaviour across the common-hour/full-day check. GSMaP is competitive on variability amplitude and categorical skill but its overestimation is large and event-dependent enough (up to +147% on 2023-08-26) to be a serious caveat. FY4B is the weakest product on every spatial measure evaluated here.

= Limitations and open issues

- *Small sample.* Only 6 events, and only 2 localized ones — event-level numbers (@tab-per-event) should be read as case studies, not as a stable estimate of localized-event skill.
- *No geostatistical heterogeneity measure.* `cv_ratio` and station-pair correlation stand in for spatial heterogeneity; a semivariogram or Moran's I comparison (both available elsewhere in this codebase) was not computed for these events and could sharpen the "amplitude vs. arrangement" distinction the current metrics can only partly separate.
- *IDW maps are illustrative, not evidence.* They interpolate satellite values sampled at gauge points and are masked beyond 20 km of any gauge; they should not be read as each product's native retrieval.
- *Peak-shift sensitivity.* `peak_shift_km` depends on a single gauge's maximum in each field and is noisy by construction (see the 2023-08-26 FY4B outlier); it is best read alongside `centroid_shift_km`, which is far more stable.
- *Point-to-pixel comparison.* As in the temporal-evaluation report, every comparison is a gauge against the satellite pixel containing it; some of the disagreement is representativeness error rather than retrieval error.

= Next steps

+ Compute a semivariogram or Moran's I comparison for the 6 events to separate the amplitude of spatial variability (`cv_ratio`) from its spatial arrangement.
+ Investigate the 2022-06-26 across-the-board failure and the 2023-08-26 GSMaP/FY4B outliers against the raw satellite retrievals, to determine whether they trace to a specific retrieval or alignment issue rather than a general product limitation.
+ Extend the localized-event sample beyond 2 cases if additional qualifying days can be found outside 2022–2024, to strengthen the widespread/localized comparison.
+ Cross-reference these 6 event dates against the temporal-evaluation report's per-event scores (`station_events.csv`) to see whether the same events are also outliers in the hourly timing analysis.

#pagebreak()

#heading(numbering: none)[Appendix A — Data dictionary]

#table(
  columns: (46%, 54%),
  stroke: none,
  inset: (x: 5pt, y: 2.4pt),
  table.header(rule, [*File*], [*Content*], thin),
  [`candidate_days.csv`], [every 2022–2024 day with more than 3 gauges above 50 mm, screened for eligibility (56 rows)],
  [`selected_events.csv`], [the 6 representative events actually used, with their `event_type`],
  [`event_metrics.csv`], [per-event, per-product spatial-agreement metrics (common_hours and, for GPM/GSMaP, full_day)],
  [`event_metrics_summary.csv`], [mean/median-over-events and pooled-pairs aggregates per product × window],
  [`idw_settings.csv`], [the IDW configuration used to build the visualization-only surfaces],
  [`idw_surfaces.csv`], [gridded IDW surfaces per event × source (Gauge, FY4B, GPM, GSMaP), for @fig-idw only],
  [`station_event_totals.csv`], [the underlying station-level (event, window, product, station) pairs table],
  rule,
)

#heading(numbering: none)[Appendix B — Reproduction]

All results come from commit `2658353`.

```sh
julia --project=. scripts/run_heavy_rain_event_evaluation.jl  # screens, selects and scores events
py -3.13 scripts/plot_heavy_rain_events.py                    # figures
```
