// Interim results: the gauge-dry station-hours - what FY4B, GPM and GSMaP, and the interpolation and fusion
// methods built on them, report when the gauges record no rain.
// Compile from this folder:  typst compile no_rain_evaluation.typ
// Every number is read from tables/norain_*.json, which scripts/plot_no_rain_evaluation.py writes from
// output/no_rain_evaluation/*.csv (scripts/run_no_rain_evaluation.jl) and output/no_rain_evaluation/fusion/*.csv
// (scripts/run_no_rain_fusion_evaluation.jl); figures are copies of output/no_rain_evaluation/figures/.

#set document(
  title: "The no-rain hours: what satellites and fused products report when the gauges are dry",
  date: datetime(year: 2026, month: 9, day: 19),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      No-rain evaluation — interim results
      #h(1fr) 19 September 2026
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
// Data. Every table is loaded once, and the prose reads single values through `val`, so a number in the
// text is always the number in the table. A lookup that finds no row stops the compile.

#let load(name) = json("tables/norain_" + name + ".json")
#let overview = load("overview")
#let anatomy = load("anatomy")
#let strata = load("strata")
#let diurnal = load("diurnal")
#let stations = load("stations")
#let baseline = load("baseline")
#let extent = load("extent")
#let spells = load("spells")
#let days = load("days")
#let daily-spells = load("dry_spells_daily")
#let calibration = load("calibration")
#let silent-spells = load("silent_spells")
#let silent-seasons = load("silent_seasons")
#let silent-sensitivity = load("silent_sensitivity")
#let screened = load("screened")
#let fusion = load("fusion")
#let fusion-strata = load("fusion_strata")
#let fusion-zero = load("fusion_zero")
#let fusion-days = load("fusion_days")
#let fusion-screened = load("fusion_screened")
#let provenance = load("provenance")

#let find(rows, ..conditions) = {
  let hits = rows.filter(row => conditions.named().pairs().all(pair => row.at(pair.at(0), default: none) == pair.at(1)))
  if hits.len() == 0 { panic("no table row for " + repr(conditions.named())) }
  hits.at(0)
}
/// One cell as typeset text: ASCII hyphens in numbers become minus signs.
#let val(rows, column, ..conditions) = find(rows, ..conditions).at(column).replace("-", "−")
/// A cell as it is, for dates and labels.
#let raw-val(rows, column, ..conditions) = find(rows, ..conditions).at(column)
#let pct(rows, column, ..conditions) = [#val(rows, column, ..conditions)%]
/// The smallest and largest value of `column` over the rows matching `keep`, as "a–b".
#let span(rows, column, keep, digits: 1) = {
  let values = rows.filter(keep).map(row => float(row.at(column).replace("+", "")))
  if values.len() == 0 { panic("empty span over " + column) }
  [#calc.round(calc.min(..values), digits: digits)–#calc.round(calc.max(..values), digits: digits)]
}

// Shorthands for the two tables the prose reads most.
#let ov(column, sample, product, threshold) = val(overview, column, sample: sample, product: product, threshold: threshold)
#let st(column, sample, product, stratifier, level) = val(strata, column, sample: sample, product: product,
  stratifier: stratifier, level: level)
#let fu(column, anchor, method, scheme: "balanced_spatial") = val(fusion, column, scheme: scheme, anchor: anchor, method: method)

/// The share of a product's false alarms (estimate >= 0.5 mm/h on a gauge-dry hour) that fall in the listed
/// levels of one stratifier: each level's share of the dry hours times its false-alarm rate.
#let fa-share(sample, product, stratifier, near) = {
  let rows = strata.filter(r => r.sample == sample and r.product == product and r.stratifier == stratifier)
  let hits = rows.filter(r => r.level in near)
  if hits.len() != near.len() { panic("fa-share: missing levels in " + stratifier) }
  let w(r) = float(r.share_of_dry) * float(r.POFD05)
  [#int(calc.round(100 * hits.map(w).sum() / rows.map(w).sum()))%]
}
/// The summed share of the dry hours (full record) in the listed levels of one stratifier.
#let share-sum(stratifier, levels) = {
  let hits = anatomy.filter(r => r.stratifier == stratifier and r.level in levels)
  if hits.len() != levels.len() { panic("share-sum: missing levels in " + stratifier) }
  calc.round(hits.map(r => float(r.share_of_dry)).sum(), digits: 1)
}
#let dd(column, sample, pair, source, threshold) = val(daily-spells, column, sample: sample, pair: pair,
  source: source, threshold: threshold)
#let sc(column, sample, product, threshold) = val(screened, column, sample: sample, product: product,
  threshold: threshold, season: "all")

#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)
#let ci(body) = text(size: 6.5pt, fill: luma(90), body.replace("-", "−"))
#let group(n, body) = table.cell(colspan: n, fill: luma(235), body)
#let note(body) = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, fill: luma(246),
  stroke: (left: 2.5pt + rgb("#2a78d6")), body)

#let P = (FY4B: "FY4B", GPM: "GPM", GSMaP: "GSMaP", GPM_uncal: "GPM uncalibrated")
#let S = (all_products: "FY4B hours", gpm_gsmap_full: "full record", gpm_cal_uncal: "full record")
#let M = (
  raw: "Raw satellite", adw: "ADW (gauges only)", idw: "IDW", tps: "TPS", gwr: "GWR (gauges only)",
  "auto": "GWR family, chosen in-fold", residual_gwr: "Residual GWR", mixed_gwr: "Mixed GWR", mgwr: "MGWR",
  blend_residual_gwr: "Residual GWR + ADW blend", blend_mixed_gwr: "Mixed GWR + ADW blend",
  blend_mgwr: "MGWR + ADW blend", blend_agrenv_residual_gwr: "Residual GWR + agreement blend",
  blend_agrenv_mixed_gwr: "Mixed GWR + agreement blend", blend_agrenv_mgwr: "MGWR + agreement blend",
  zero: "Always zero", train_clim: "Training climatology", hour_field_mean: "Hourly field mean",
)
#let A = (FY4B: "FY4B", GPM: "GPM", GSMaP: "GSMaP", MERGED_MEAN: "Merged (equal)", MERGED_OLS: "Merged (OLS)")

// ---------------------------------------------------------------------------------------------
// Title

#align(left)[
  #set text(hyphenate: false)
  #text(font: ("Segoe UI", "Arial"), size: 9pt, fill: luma(100), tracking: 0.04em)[INTERIM RESULTS · 19 SEPTEMBER 2026]
  #v(0.2em)
  #text(font: ("Segoe UI", "Arial"), size: 19pt, weight: "semibold")[What do the products report when it does not rain?]
  #v(0.1em)
  #text(size: 12pt, fill: luma(60))[The #val(anatomy, "dry_share", stratifier: "all", level: "all")% of gauge station-hours with no rain: FY4B, GPM IMERG (calibrated and uncalibrated), GSMaP and the gauge–satellite fusion methods at 237 hourly gauges, 2022–2024]
]
#v(0.6em)
#line(length: 100%, stroke: 0.5pt + luma(180))

= Key findings

- *The products rain on 6–7% of dry hours, mostly as drizzle.* Over 2022–2024 GPM reports ≥ 0.1 mm/h on #ov("POFD", "gpm_gsmap_full", "GPM", "0.1")% of gauge-dry hours and GSMaP on #ov("POFD", "gpm_gsmap_full", "GSMaP", "0.1")%, with too many wet hours at 0.1 mm/h (frequency bias #ov("freq_bias", "gpm_gsmap_full", "GPM", "0.1"), #ov("freq_bias", "gpm_gsmap_full", "GSMaP", "0.1")) but too few at the gauge's 0.5 mm resolution (#ov("freq_bias", "gpm_gsmap_full", "GPM", "0.5"), #ov("freq_bias", "gpm_gsmap_full", "GSMaP", "0.5")). A third of their rain volume falls on gauge-dry hours; for FY4B it is #ov("dry_volume_share", "all_products", "FY4B", "0.1")%.
- *Most of that rain is displaced, not invented.* GPM's false-alarm rate at ≥ 0.5 mm/h is #st("POFD05", "gpm_gsmap_full", "GPM", "rain_proximity", "1 h")% within an hour of the gauge's own rain and #st("POFD05", "gpm_gsmap_full", "GPM", "rain_proximity", "> 24 h")% more than a day from it; #fa-share("gpm_gsmap_full", "GPM", "rain_proximity", ("1 h", "2-3 h", "4-6 h")) of its false alarms are within 6 h of rain. At 0.5 mm/h GPM and GSMaP false-alarm about as often as a neighbouring gauge 15–30 km away does. FY4B is the exception: its false alarms are the least tied to real rain, and it produces whole phantom rain areas.
- *The products share their false alarms,* so merging them cannot remove them; they rain early, before the gauges' onset; and in dry air most of their wet hours are false (false-alarm ratio #st("FAR01", "gpm_gsmap_full", "GPM", "dewpoint_depression", ">= 10")% at a dew-point depression ≥ 10 °C), the signature of rain evaporating below the cloud.
- *IMERG's gauge calibration adds dry-hour rain rather than removing it* (spurious rain #ov("spurious_mm", "gpm_cal_uncal", "GPM_uncal", "0.1") → #ov("spurious_mm", "gpm_cal_uncal", "GPM", "0.1") mm per year): it rescales amounts and leaves the occurrence pattern as it was.
- *Some gauges are silent.* #val(silent-sensitivity, "n_spells", min_wet_days: "3") spells at #val(silent-sensitivity, "n_stations", min_wet_days: "3") gauges read 0.0 for a week to more than six months while their neighbours recorded rain. They are #val(silent-sensitivity, "dry_hour_share", min_wet_days: "3")% of the dry hours and move no pooled score, but the benchmark trains and scores on their zeros.
- *For the fusion methods dry hours decide the ranking.* Gauge-only ADW is no drier than GPM on held-out dry hours (≥ 0.1 mm/h on #fu("POFD01", "GPM", "adw")% against #fu("POFD01", "GPM", "raw")%), but its errors are small. Unblended MGWR copies the anchor's false rain and its dry-hour RMSE is #fu("dry_vs_adw", "GPM", "mgwr").replace("−", "")% worse than ADW's; the blend toward ADW repairs this and beats ADW on GPM's dry hours by #fu("dry_vs_adw", "GPM", "blend_mgwr")%. The OLS-merged anchor never predicts a dry day, yet has the lowest dry-hour RMSE.

= Question and scope

The three earlier interim reports looked at rain: its intensity and timing (temporal evaluation), its spatial pattern on heavy-rain days (heavy-rain events), and what gauges and satellites gain from each other (fusion). Yet #val(anatomy, "dry_share", stratifier: "all", level: "all")% of gauge station-hours record no rain at all, and the fusion report found that dry hours decide how the methods rank: unblended fusion loses to gauge-only interpolation mainly on hours when the gauge is dry but the satellite is not, and the blending rule wins overall mainly by repairing those hours. This report looks at the no-rain hours on their own, as widely as the data allow:

+ How often, how much and when do the products report rain when the gauge reports none (@sec-how-often, @sec-when-where)?
+ Is that rain really false — or is it rain next to the gauge in time or space, rain below the gauge's resolution, or rain the gauge failed to record (@sec-false, @sec-silent)?
+ Do the products reproduce dry spells and dry days (@sec-spells)? Does IMERG's gauge calibration change any of this (@sec-calibration)?
+ What do the interpolation and fusion methods predict on dry hours, and where does their dry-hour error come from (@sec-fusion)?

= Data and samples

- *Gauges.* 237 hourly tipping-bucket gauges (0.5 mm resolution) in the Shiyan study area, hour-ending Beijing time (BJT), 08–08 BJT days from 1 January 2022 to 31 December 2024.
- *Products*, at the gauge pixel and on the gauge clock: FY4B (strict hourly aggregation, navigation-corrected), GPM IMERG V07 Final `precipitation` (gauge-calibrated) and GSMaP V8 operational (no gauge correction), as in the earlier reports. New here: IMERG V07 `precipitationUncal`, the same multi-satellite estimate *before* its monthly gauge calibration, downloaded half-hourly from NASA GES DISC and placed on the gauge clock (@sec-calibration).
- *ERA5-Land* 2 m temperature and dew point at every gauge, hourly, on the same clock the benchmark's covariates use.
- *Three samples*, as in the temporal report (every product in a sample is scored on identical station-hours):
  - *FY4B hours*: cells where the gauge and all three products report; FY4B never has the 23:00 label and misses whole months.
  - *Full record*: GPM and GSMaP over all 26,304 hours.
  - *Benchmark grid*: the 13,471 hours of the interpolation benchmark, for the fusion methods (@sec-fusion), each anchor on its own evaluation mask.
- *Fusion methods.* The held-out predictions the full benchmark run saved for every method, anchor and cross-validation scheme (`oof_<method>.csv`): nothing is refitted. The run is the one the fusion report used.

= Methods

== What counts as no rain

A station-hour is *gauge-dry* when the gauge reads below 0.1 mm — for a 0.5 mm bucket, a reading of 0.0.#footnote[A handful of gauge readings are non-zero but below 0.1 mm; they are counted as dry, exactly as the benchmark counts them.] A product is *wet* at a threshold $t$ when it reaches $t$; two thresholds are used throughout, 0.1 mm/h (the wet threshold) and 0.5 mm/h (the gauge resolution), because a product value between them may be drizzle the bucket had not yet tipped for. With the 2×2 table of gauge and product occurrence (@tab-scores), the dry side is described by:

#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    stroke: none,
    inset: (x: 6pt, y: 2.6pt),
    table.header(rule, [*Score*], [*Meaning*], thin),
    [False-alarm rate (POFD)], [P(product wet | gauge dry) — how often a dry gauge hour is called wet; $1 − "POFD"$ is the specificity],
    [False-alarm ratio (FAR)], [P(gauge dry | product wet) — what share of the product's wet hours are dry at the gauge],
    [Dry reliability (NPV)], [P(gauge dry | product dry) — when the product says dry, how often the gauge agrees],
    [Frequency bias], [product wet hours / gauge wet hours; above 1 the product rains too often],
    [PSS, HSS], [Peirce and Heidke skill scores of the whole 2×2 table (1 perfect, 0 no skill)],
    [Dry-hour amount], [the product's mean on gauge-dry hours (mm/h); its share of the product's total volume; the share of gauge-dry hours on which it is exactly 0],
    [Spurious rain], [the product's total on gauge-dry hours per year of scored record (mm/yr). Only meaningful for the full record: the FY4B-hours sample is weighted to summer],
    rule,
  ),
  caption: [Scores of the no-rain side of the occurrence table.],
) <tab-scores>
]

Uncertainty is a *day-block bootstrap* (1,000 replicates, fixed seed): whole 08–08 days are resampled, keeping the hours and gauges of one weather system together.

== Where a dry hour sits

A satellite value on a dry gauge hour can be wrong in four quite different ways, and the report separates them with four descriptions of the gauge network around each dry hour:

- *Distance to the gauge's own rain*: hours to the nearest wet hour at the same gauge, before or after (1 h, 2–3 h, 4–6 h, 7–24 h, > 24 h). Rain an hour away is a timing error; rain a day away is not.
- *Side of the rain* (± 3 h): a gap inside an event (rain on both sides), just before the onset, just after the end, or away from rain.
- *Neighbouring gauges* within 15 km: all dry, fewer than half wet, or half or more wet. A satellite pixel is an area; rain at a neighbour is rain inside or near the pixel.
- *The whole network*: the share of reporting gauges that are wet in that hour — network-dry (none), isolated (≤ 5%), scattered (5–25%) or widespread (> 25%). A network-dry hour is as close to "no rain anywhere" as the gauges can show.

Hours whose class depends on a missing gauge value are left out of that stratifier.

== Two physical baselines

- *A gauge as the estimate.* Every pair of gauges within 30 km is scored as if one were a product estimating the other: how often does gauge B read ≥ 0.1 mm while gauge A is dry, by separation? Two gauges a few km apart sample one rain field at two points, so this is the false-alarm rate that point-to-area comparison produces even for a perfect estimate — the floor under the satellites' rates.
- *The air at the gauge.* ERA5-Land dew-point depression $T − T_d$ (how far the near-surface air is from saturation) and 2 m temperature. Rain that evaporates below the cloud base, and retrievals over cold or snow-covered surfaces, are the two textbook sources of satellite rain that never reaches a gauge.

== The silent-gauge screen <sec-silent-method>

A dry hour is only as good as the gauge. A gauge that reads 0.0 while it is not measuring — a blocked funnel, a stopped logger — produces dry hours on which any rain at all counts as a satellite false alarm. The screen works on 08–08 daily totals: a run of at least 7 complete days below 1 mm at a gauge is flagged when the median daily total of its three nearest reporting gauges reaches 5 mm on at least 3 of those days. The rule rests on one measured rate: over all station-days, a gauge stays below 1 mm on only #val(silent-sensitivity, "p_dry_given_wet", min_wet_days: "3")% of the days its neighbours' median reaches 5 mm, so three such days inside one dry run are very unlikely to be weather. The sensitivity to that "3" is reported (@tab-silent-rule). The screen is a description, not a correction: no gauge value is changed, and every result elsewhere is on the unscreened data unless stated.

== Dry spells and dry days

An *hourly dry spell* is a run of consecutive hours below 0.1 mm (products also at 0.5 mm/h). Spells cut by a missing hour are censored and not counted, so hourly spells are only computed on the full record: FY4B's daily 23:00 gap would cut every spell at under a day. *Dry days* use 08–08 BJT totals at 0.1 mm (the Chinese "no-rain day") and 1 mm (the ETCCDI dry day). On the FY4B hours and the benchmark grid a day counts when at least 20 of its 24 hours are scored, and gauge and product are summed over exactly the same hours.

== Interpolation and fusion methods

The fusion methods are those of the fusion report: gauge-only interpolation (IDW, ADW, TPS, GWR on coordinates), residual correction of an anchor (residual, mixed and multiscale GWR, and the in-fold choice among them), and the two blends toward ADW (constant, and by product agreement), on five anchors (FY4B, GPM, GSMaP and their equal and OLS merges). Each method is scored on the held-out gauges' dry hours; balanced spatial cross-validation is the primary scheme, random CV the robustness check. Paired day-block bootstraps compare each method with ADW on identical cells, and the squared-error gap to ADW is split into four quadrants: gauge dry or wet × anchor dry or wet. Before any table is written, every method's dry-hour n, RMSE and bias are checked against the benchmark run's own pooled metrics.

= Results

== The no-rain record <sec-record>

Over the full record #val(anatomy, "dry_share", stratifier: "all", level: "all")% of station-hours are gauge-dry (@fig-anatomy). The share hardly changes between spring, summer and autumn (#val(anatomy, "dry_share", stratifier: "season", level: "SON")–#val(anatomy, "dry_share", stratifier: "season", level: "JJA")%) and rises to #val(anatomy, "dry_share", stratifier: "season", level: "DJF")% in winter; by month it runs from #val(anatomy, "dry_share", stratifier: "month", level: "Oct")% in October to #val(anatomy, "dry_share", stratifier: "month", level: "Dec")% in December. Over the day the driest hour is #raw-val(anatomy, "share_of_dry", stratifier: "hour_extremes", level: "driest") BJT (#val(anatomy, "dry_share", stratifier: "hour_extremes", level: "driest")% dry) and the wettest #raw-val(anatomy, "share_of_dry", stratifier: "hour_extremes", level: "wettest") BJT (#val(anatomy, "dry_share", stratifier: "hour_extremes", level: "wettest")%). The relief regions differ by little more than two points.

Most dry hours are far from any rain. #val(anatomy, "share_of_dry", stratifier: "rain_proximity", level: "> 24 h")% of them are more than a day from the gauge's own nearest wet hour and only #val(anatomy, "share_of_dry", stratifier: "rain_proximity", level: "1 h")% are within an hour of it; in #val(anatomy, "share_of_dry", stratifier: "neighbours", level: "all neighbours dry")% every gauge within 15 km is dry as well, and #val(anatomy, "share_of_dry", stratifier: "network", level: "network dry")% fall in hours when no reporting gauge in the whole network is wet. Such network-dry hours are #val(anatomy, "dry_share", stratifier: "network_hours", level: "network dry")% of the #raw-val(anatomy, "gauge_mm", stratifier: "network_hours", level: "network dry") hours with at least 90% of the gauges reporting. By count, then, a satellite's dry-hour record is dominated by clear weather; the next sections show that its false alarms are not.

#figure(
  image("figures/norain_fig1_anatomy.png", width: 100%),
  caption: [The no-rain record at 237 gauges, full 2022–2024 record: dry share by month and by hour and season, and how the dry hours split by distance to the gauge's own rain and by the state of the whole network.],
) <fig-anatomy>

== How often the products rain on dry hours <sec-how-often>

@tab-overview lists the scores and @fig-exceedance the rate at every threshold. Over the full record GPM reports ≥ 0.1 mm/h on #ov("POFD", "gpm_gsmap_full", "GPM", "0.1")% of gauge-dry hours and GSMaP on #ov("POFD", "gpm_gsmap_full", "GSMaP", "0.1")%; at the gauge's own resolution, ≥ 0.5 mm/h, the rates fall to #ov("POFD", "gpm_gsmap_full", "GPM", "0.5")% and #ov("POFD", "gpm_gsmap_full", "GSMaP", "0.5")%. The frequency bias shows what kind of rain this is: at 0.1 mm/h both products have too many wet hours (#ov("freq_bias", "gpm_gsmap_full", "GPM", "0.1") and #ov("freq_bias", "gpm_gsmap_full", "GSMaP", "0.1")), at 0.5 mm/h too few (#ov("freq_bias", "gpm_gsmap_full", "GPM", "0.5") and #ov("freq_bias", "gpm_gsmap_full", "GSMaP", "0.5")). The surplus of wet hours is drizzle between 0.1 and 0.5 mm/h, which a 0.5 mm bucket cannot record in one hour. More than half of the products' wet hours at 0.1 mm/h are gauge-dry (FAR #ov("FAR", "gpm_gsmap_full", "GPM", "0.1")% and #ov("FAR", "gpm_gsmap_full", "GSMaP", "0.1")%).

The rain is small per hour but large in total. On gauge-dry hours GPM and GSMaP are exactly zero #ov("zero_share", "gpm_gsmap_full", "GPM", "0.1")% and #ov("zero_share", "gpm_gsmap_full", "GSMaP", "0.1")% of the time, yet #ov("dry_volume_share", "gpm_gsmap_full", "GPM", "0.1")% of GPM's and #ov("dry_volume_share", "gpm_gsmap_full", "GSMaP", "0.1")% of GSMaP's rain falls on those hours: #ov("spurious_mm", "gpm_gsmap_full", "GPM", "0.1") and #ov("spurious_mm", "gpm_gsmap_full", "GSMaP", "0.1") mm per year, against a gauge total of #ov("gauge_mm", "gpm_gsmap_full", "GPM", "0.1") mm per year.

FY4B behaves differently. On the hours it covers it reports ≥ 0.1 mm/h on fewer dry hours than GPM and GSMaP (#ov("POFD", "all_products", "FY4B", "0.1")% against #ov("POFD", "all_products", "GPM", "0.1")% and #ov("POFD", "all_products", "GSMaP", "0.1")%) but ≥ 0.5 mm/h on more (#ov("POFD", "all_products", "FY4B", "0.5")% against #ov("POFD", "all_products", "GPM", "0.5")% and #ov("POFD", "all_products", "GSMaP", "0.5")%), and its exceedance curve is flat below 0.1 mm/h (@fig-exceedance): FY4B has almost no light values, so when it is wet on a dry hour it is wet with a real amount. Its mean on dry hours (#ov("mean_dry", "all_products", "FY4B", "0.1") mm/h) is twice GPM's (#ov("mean_dry", "all_products", "GPM", "0.1")), and #ov("dry_volume_share", "all_products", "FY4B", "0.1")% of its rain volume falls on gauge-dry hours, against #ov("dry_volume_share", "all_products", "GPM", "0.1")% for GPM and #ov("dry_volume_share", "all_products", "GSMaP", "0.1")% for GSMaP.

#figure(
  image("figures/norain_fig2_exceedance.png", width: 100%),
  caption: [Share of gauge-dry station-hours on which each product reaches the threshold $t$. Dotted: a neighbouring gauge 0–5 km away used as the estimate.],
) <fig-exceedance>

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.4pt),
    table.header(rule, [*Sample*], [*Product*], [*$t$*], [*POFD %*], [*FAR %*], [*NPV %*], [*Freq. bias*], [*HSS*],
      [*Dry mean*], [*Dry vol. %*], thin),
    ..overview.map(r => (
      S.at(r.sample), P.at(r.product), r.threshold,
      [#r.POFD #ci(r.POFD_ci)], r.FAR, r.NPV, [#r.freq_bias #ci(r.freq_bias_ci)], r.HSS, r.mean_dry,
      r.dry_volume_share,
    )).flatten(),
    rule,
  ),
  caption: [Occurrence scores on gauge-dry station-hours; gauge wet at ≥ 0.1 mm/h, estimate wet at $t$ (mm/h). Brackets: 95% day-block bootstrap. Dry mean: the estimate's mean on gauge-dry hours (mm/h). Dry vol.: share of the estimate's total falling on gauge-dry hours. The full-record GPM rows appear twice because the calibrated/uncalibrated sample is scored separately.],
) <tab-overview>
]

The gauges' own spread sets the scale for these numbers (@fig-floor). Scored against one another, a gauge 0–5 km away reports rain on #val(baseline, "POFD", sample: "gpm_gsmap_full", level: "0-5 km")% of the hours its neighbour is dry, and one 20–30 km away on #val(baseline, "POFD", sample: "gpm_gsmap_full", level: "20-30 km")%. At ≥ 0.5 mm/h GPM and GSMaP (#ov("POFD", "gpm_gsmap_full", "GPM", "0.5")% and #ov("POFD", "gpm_gsmap_full", "GSMaP", "0.5")%) sit inside that range: at the gauge resolution they false-alarm about as often as a gauge 15–30 km away would. At 0.1 mm/h they do so two to four times as often. FY4B (#ov("POFD", "all_products", "FY4B", "0.5")% on its hours) is above every gauge baseline.

Across gauges the rates are uniform. At ≥ 0.5 mm/h, 80% of the gauges lie between #val(stations, "POFD_p10", sample: "gpm_gsmap_full", product: "GPM", threshold: "0.5")% and #val(stations, "POFD_p90", sample: "gpm_gsmap_full", product: "GPM", threshold: "0.5")% for GPM and between #val(stations, "POFD_p10", sample: "all_products", product: "FY4B", threshold: "0.5")% and #val(stations, "POFD_p90", sample: "all_products", product: "FY4B", threshold: "0.5")% for FY4B (@fig-maps), and the three relief regions are within a few tenths of a point of each other. The false alarms are not the problem of a few pixels or of the mountains.

#figure(
  image("figures/norain_fig3_contingency.png", width: 100%),
  caption: [The dry side of the occurrence table for every product and sample, at 0.1 and 0.5 mm/h, with 95% day-block bootstrap intervals.],
) <fig-contingency>

#figure(
  image("figures/norain_fig5_station_maps.png", width: 88%),
  caption: [False-alarm rate at ≥ 0.5 mm/h (top) and the share of the product's rain that falls on gauge-dry hours (bottom), at each gauge, on the hours FY4B covers.],
) <fig-maps>

== When and where <sec-when-where>

*Season.* The false alarms follow the rain season. GPM's rate at ≥ 0.5 mm/h is #st("POFD05", "gpm_gsmap_full", "GPM", "season", "JJA")% in summer and #st("POFD05", "gpm_gsmap_full", "GPM", "season", "DJF")% in winter. In winter the few wet hours the products do report are mostly gauge-dry: at 0.1 mm/h the false-alarm ratio is #st("FAR01", "gpm_gsmap_full", "GPM", "season", "DJF")% for GPM and #st("FAR01", "gpm_gsmap_full", "GSMaP", "season", "DJF")% for GSMaP, against #st("FAR01", "gpm_gsmap_full", "GPM", "season", "JJA")% and #st("FAR01", "gpm_gsmap_full", "GSMaP", "season", "JJA")% in summer. Part of this is the gauge rather than the satellite: a tipping bucket records snow only when it melts and nothing while frozen, and winter is where the silent-gauge screen finds many of its spells (@sec-silent). FY4B hardly rains in winter at all — #st("POFD01", "all_products", "FY4B", "season", "DJF")% of dry hours at 0.1 mm/h, with a winter frequency bias of #st("freq_bias01", "all_products", "FY4B", "season", "DJF") — and below 0 °C essentially never (#st("POFD01", "all_products", "FY4B", "t2m", "< 0")%): it has few false alarms in the cold season because it misses the cold-season rain as well.

*Hour of day* (@fig-diurnal). GPM's false alarms peak at #val(diurnal, "peak_hour", sample: "gpm_gsmap_full", product: "GPM", season: "all") BJT (#val(diurnal, "peak_POFD05", sample: "gpm_gsmap_full", product: "GPM", season: "all")% at ≥ 0.5 mm/h; #val(diurnal, "peak_POFD05", sample: "gpm_gsmap_full", product: "GPM", season: "JJA")% in summer), with the afternoon convection the gauges record, and are lowest at #val(diurnal, "low_hour", sample: "gpm_gsmap_full", product: "GPM", season: "all") BJT. GSMaP's peak at #val(diurnal, "peak_hour", sample: "gpm_gsmap_full", product: "GSMaP", season: "all") BJT (#val(diurnal, "peak_hour", sample: "gpm_gsmap_full", product: "GSMaP", season: "MAM") BJT in spring, #val(diurnal, "peak_hour", sample: "gpm_gsmap_full", product: "GSMaP", season: "JJA") in summer) has no counterpart in the gauges; it is the same spurious late-morning maximum the temporal report found in GSMaP's rain. FY4B's false alarms peak in the evening (#val(diurnal, "peak_hour", sample: "all_products", product: "FY4B", season: "all") BJT on its hours).

#figure(
  image("figures/norain_fig4_diurnal.png", width: 100%),
  caption: [False-alarm rate at ≥ 0.5 mm/h by hour of day: all seasons on the hours FY4B covers (left), and by season over the full record.],
) <fig-diurnal>

== Is the rain really false? <sec-false>

A gauge-dry hour is a statement about one 200 cm² funnel in one hour. @fig-context sorts the dry hours by what the gauges show around them, and the false-alarm rate changes by more than an order of magnitude between the groups (full record, GPM, ≥ 0.5 mm/h):

- *In time.* #st("POFD05", "gpm_gsmap_full", "GPM", "rain_proximity", "1 h")% within an hour of the gauge's own rain, #st("POFD05", "gpm_gsmap_full", "GPM", "rain_proximity", "> 24 h")% more than a day from it. Weighted by how many dry hours each group holds, #fa-share("gpm_gsmap_full", "GPM", "rain_proximity", ("1 h", "2-3 h", "4-6 h")) of GPM's false alarms fall within 6 h of rain at the same gauge (GSMaP #fa-share("gpm_gsmap_full", "GSMaP", "rain_proximity", ("1 h", "2-3 h", "4-6 h")), FY4B on its hours #fa-share("all_products", "FY4B", "rain_proximity", ("1 h", "2-3 h", "4-6 h"))), although those hours are only #share-sum("rain_proximity", ("1 h", "2-3 h", "4-6 h"))% of the dry record.
- *Before or after.* #st("POFD05", "gpm_gsmap_full", "GPM", "rain_side", "before rain")% in the three hours before the rain starts against #st("POFD05", "gpm_gsmap_full", "GPM", "rain_side", "after rain")% in the three hours after it ends (GSMaP #st("POFD05", "gpm_gsmap_full", "GSMaP", "rain_side", "before rain")% and #st("POFD05", "gpm_gsmap_full", "GSMaP", "rain_side", "after rain")%): the satellites start the rain early, as the temporal report found for its onset. Dry gaps inside an event are called wet #st("POFD05", "gpm_gsmap_full", "GPM", "rain_side", "inside an event gap")% of the time.
- *In space.* #st("POFD05", "gpm_gsmap_full", "GPM", "neighbours", "half or more wet")% when at least half the gauges within 15 km are wet, #st("POFD05", "gpm_gsmap_full", "GPM", "neighbours", "all neighbours dry")% when all are dry. On network-dry hours the rate is #st("POFD05", "gpm_gsmap_full", "GPM", "network", "network dry")%, on hours when more than a quarter of the network is wet #st("POFD05", "gpm_gsmap_full", "GPM", "network", "widespread (> 25% wet)")%.

Most of GPM's and GSMaP's "false" rain is therefore real rain in the wrong place or at the wrong time by a few kilometres or hours — the point-to-pixel and timing error of a correct rain field — not rain from nothing.

#figure(
  image("figures/norain_fig6_context.png", width: 100%),
  caption: [False-alarm rate at ≥ 0.5 mm/h on gauge-dry hours grouped by the gauges around them, on the hours FY4B covers. In brackets: each group's share of the dry hours.],
) <fig-context>

*FY4B's false alarms are the least tied to real rain.* On the hours it covers FY4B is wet on #st("POFD05", "all_products", "FY4B", "network", "network dry")% of network-dry hours, against #st("POFD05", "all_products", "GPM", "network", "network dry")% for GPM and #st("POFD05", "all_products", "GSMaP", "network", "network dry")% for GSMaP, and #fa-share("all_products", "FY4B", "network", ("network dry",)) of its false alarms fall on such hours (GPM #fa-share("all_products", "GPM", "network", ("network dry",))). On network-dry hours FY4B covers more than half of the gauge pixels on #val(extent, "gt50", sample: "all_products", product: "FY4B")% of the hours, GPM on #val(extent, "gt50", sample: "all_products", product: "GPM")% (@fig-agreement, right): FY4B produces whole phantom rain areas, the passive-microwave products scattered pixels. GPM reports rain at one gauge pixel or more on #val(extent, "any", sample: "all_products", product: "GPM")% of network-dry hours but at more than 5% of them on only #val(extent, "gt5", sample: "all_products", product: "GPM")%.

*The products' false alarms are shared* (@fig-agreement). When neither of the other two products is wet, a product is wet at ≥ 0.5 mm/h on under 2% of gauge-dry hours (GPM #st("POFD05", "all_products", "GPM", "other_products", "0 others wet")%); when both others are wet, GPM is wet on #st("POFD05", "all_products", "GPM", "other_products", "2 others wet")% and GSMaP on #st("POFD05", "all_products", "GSMaP", "other_products", "2 others wet")%. Three products built from different sensors agreeing on rain the gauge did not record is further evidence that much of it is rain near the gauge rather than retrieval noise — and it means that averaging the products cannot remove it.

#figure(
  image("figures/norain_fig8_agreement.png", width: 100%),
  caption: [Product agreement on gauge-dry hours (left, middle), and how widely each product rains on hours when every reporting gauge is dry (right). Hours FY4B covers.],
) <fig-agreement>

*The air at the gauge* (@fig-floor, middle and right). The false-alarm rate falls as the near-surface air gets drier: GPM's rate at ≥ 0.5 mm/h is #st("POFD05", "gpm_gsmap_full", "GPM", "dewpoint_depression", "< 1")% at a dew-point depression below 1 °C and #st("POFD05", "gpm_gsmap_full", "GPM", "dewpoint_depression", ">= 10")% at 10 °C or more. The gauge's rain falls much faster, so in dry air the few wet hours the products report are mostly false: at a depression of 10 °C or more the false-alarm ratio at 0.1 mm/h is #st("FAR01", "gpm_gsmap_full", "GPM", "dewpoint_depression", ">= 10")% for GPM and #st("FAR01", "gpm_gsmap_full", "GSMaP", "dewpoint_depression", ">= 10")% for GSMaP, with frequency biases of #st("freq_bias01", "gpm_gsmap_full", "GPM", "dewpoint_depression", ">= 10") and #st("freq_bias01", "gpm_gsmap_full", "GSMaP", "dewpoint_depression", ">= 10")\; near saturation the bias is #st("freq_bias01", "gpm_gsmap_full", "GPM", "dewpoint_depression", "< 1"). This is the signature of rain that evaporates below the cloud base: the satellite sees precipitation aloft that never reaches the ground. Warm air shows the same (above 25 °C: false-alarm ratio #st("FAR01", "gpm_gsmap_full", "GPM", "t2m", ">= 25")%, bias #st("freq_bias01", "gpm_gsmap_full", "GPM", "t2m", ">= 25")).

#figure(
  image("figures/norain_fig7_floor_meteorology.png", width: 100%),
  caption: [Left: how often a gauge reports rain when a neighbour at the given separation is dry, against each product's rate at ≥ 0.5 mm/h. Middle, right: the products' false-alarm rate by ERA5-Land dew-point depression and 2 m temperature at the gauge. Hours FY4B covers.],
) <fig-floor>

== Silent gauges <sec-silent>

Some gauges read 0.0 for weeks or months while the gauges around them record rain. The screen of @sec-silent-method flags #val(silent-sensitivity, "n_spells", min_wet_days: "3") such spells at #val(silent-sensitivity, "n_stations", min_wet_days: "3") of the 237 gauges: #val(silent-sensitivity, "station_days", min_wet_days: "3") station-days, holding #val(silent-sensitivity, "dry_hours", min_wet_days: "3") gauge-dry hours (#val(silent-sensitivity, "dry_hour_share", min_wet_days: "3")% of all of them). The longest (@tab-silent-spells) are unmistakable: gauge #silent-spells.at(0).station_id logged #silent-spells.at(0).gauge_total mm over #silent-spells.at(0).n_days days from #silent-spells.at(0).first_date to #silent-spells.at(0).last_date, while the median of its three nearest gauges summed to #silent-spells.at(0).neighbour_total mm and GPM to #silent-spells.at(0).GPM_total mm. Such a spell is also the gauge record's longest hourly dry spell, #val(spells, "max_length", sample: "gpm_gsmap_full", source: "Gauge", threshold: "0.1") h, and the gauges at the far right of @fig-spells (right). The flagged days fall in every season (spring #val(silent-seasons, "share", season: "MAM")%, summer #val(silent-seasons, "share", season: "JJA")%, autumn #val(silent-seasons, "share", season: "SON")%, winter #val(silent-seasons, "share", season: "DJF")%), so frozen buckets explain only part of them.

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.4pt),
    table.header(rule, [*Gauge*], [*Spell*], [*Days*], [*Gauge mm*], [*Neighbours mm*], [*GPM mm*], [*GSMaP mm*], thin),
    ..silent-spells.slice(0, 8).map(r => (
      r.station_id, [#r.first_date – #r.last_date], r.n_days, r.gauge_total, r.neighbour_total, r.GPM_total, r.GSMaP_total,
    )).flatten(),
    rule,
  ),
  caption: [The eight flagged spells with the most rain at the neighbouring gauges. Neighbours: the sum over the spell of the daily median of the three nearest reporting gauges.],
) <tab-silent-spells>
]

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right),
    stroke: none,
    inset: (x: 6pt, y: 2.4pt),
    table.header(rule, [*Neighbour-wet days required*], [*Spells*], [*Gauges*], [*Station-days*], [*Dry hours*], [*Share of dry hours*], thin),
    ..silent-sensitivity.map(r => (r.min_wet_days, r.n_spells, r.n_stations, r.station_days, r.dry_hours, [#r.dry_hour_share%])).flatten(),
    rule,
  ),
  caption: [Sensitivity of the silent-gauge screen to the number of days, within a dry run of at least 7 days, on which the neighbours' median reaches 5 mm. A working gauge stays below 1 mm on #val(silent-sensitivity, "p_dry_given_wet", min_wet_days: "3")% of such days.],
) <tab-silent-rule>
]

Inside the flagged spells every product "false-alarms" far more often (@fig-silent): over the full record GPM at ≥ 0.5 mm/h #sc("POFD_silent", "gpm_gsmap_full", "GPM", "0.5")% #ci(sc("POFD_silent_ci", "gpm_gsmap_full", "GPM", "0.5")) against #sc("POFD_working", "gpm_gsmap_full", "GPM", "0.5")% #ci(sc("POFD_working_ci", "gpm_gsmap_full", "GPM", "0.5")) at working gauges, GSMaP #sc("POFD_silent", "gpm_gsmap_full", "GSMaP", "0.5")% against #sc("POFD_working", "gpm_gsmap_full", "GSMaP", "0.5")%. FY4B shows the smallest difference (#sc("POFD_silent", "all_products", "FY4B", "0.5")% against #sc("POFD_working", "all_products", "FY4B", "0.5")%), as expected if its false alarms are the least tied to real rain. But the spells are too few to move the pooled scores: they carry #sc("silent_false_alarm_share", "gpm_gsmap_full", "GPM", "0.5")% of GPM's false alarms, and removing them changes its false-alarm rate from #sc("POFD_all", "gpm_gsmap_full", "GPM", "0.5")% to #sc("POFD_working", "gpm_gsmap_full", "GPM", "0.5")% and its spurious rain from #sc("spurious_all", "gpm_gsmap_full", "GPM", "0.1") to #sc("spurious_working", "gpm_gsmap_full", "GPM", "0.1") mm per year. Silent gauges are a data-quality problem worth fixing — the benchmark trains on their zeros — but they do not explain the satellites' dry-hour rain.

#figure(
  image("figures/norain_fig13_silent_gauges.png", width: 100%),
  caption: [Left: the flagged silent-gauge spells in time, one row per gauge. Right: the products' false-alarm rate at ≥ 0.5 mm/h on dry hours at working gauges and inside the flagged spells.],
) <fig-silent>

== Dry spells and dry days <sec-spells>

*Hourly spells* (full record, @fig-spells left). At 0.1 mm/h the products break the gauges' long dry spells with drizzle: the mean complete dry spell is #val(spells, "mean_length", sample: "gpm_gsmap_full", source: "Gauge", threshold: "0.1") h at the gauges, #val(spells, "mean_length", sample: "gpm_gsmap_full", source: "GPM", threshold: "0.1") h for GPM and #val(spells, "mean_length", sample: "gpm_gsmap_full", source: "GSMaP", threshold: "0.1") h for GSMaP, and the share of dry hours inside spells of three days or more is #val(spells, "share_in_ge72", sample: "gpm_gsmap_full", source: "Gauge", threshold: "0.1")% against #val(spells, "share_in_ge72", sample: "gpm_gsmap_full", source: "GPM", threshold: "0.1")% and #val(spells, "share_in_ge72", sample: "gpm_gsmap_full", source: "GSMaP", threshold: "0.1")%. At 0.5 mm/h the products' spells are as long as the gauges' or longer (#val(spells, "mean_length", sample: "gpm_gsmap_full", source: "GPM", threshold: "0.5") and #val(spells, "mean_length", sample: "gpm_gsmap_full", source: "GSMaP", threshold: "0.5") h; #val(spells, "share_in_ge72", sample: "gpm_gsmap_full", source: "GPM", threshold: "0.5")% in spells of three days or more) — the same drizzle story as the frequency bias.

*Dry days.* Over the full record #dd("dry_share", "gpm_gsmap_full", "GPM", "Gauge", "0.1")% of gauge days are below 0.1 mm, against #dd("dry_share", "gpm_gsmap_full", "GPM", "GPM", "0.1")% of GPM days and #dd("dry_share", "gpm_gsmap_full", "GSMaP", "GSMaP", "0.1")% of GSMaP days; the wet-day frequency bias is #val(days, "freq_bias", sample: "gpm_gsmap_full", product: "GPM", stratifier: "all", level: "all", threshold: "0.1") and #val(days, "freq_bias", sample: "gpm_gsmap_full", product: "GSMaP", stratifier: "all", level: "all", threshold: "0.1"). At the ETCCDI 1 mm definition the gap narrows (#dd("dry_share", "gpm_gsmap_full", "GPM", "Gauge", "1")% of gauge days dry, #dd("dry_share", "gpm_gsmap_full", "GPM", "GPM", "1")% for GPM; bias #val(days, "freq_bias", sample: "gpm_gsmap_full", product: "GPM", stratifier: "all", level: "all", threshold: "1")). On the hours it covers, FY4B is the closest to the gauges in wet-day frequency (#val(days, "freq_bias", sample: "all_products", product: "FY4B", stratifier: "all", level: "all", threshold: "0.1") at 0.1 mm, #val(days, "freq_bias", sample: "all_products", product: "FY4B", stratifier: "all", level: "all", threshold: "1") at 1 mm, against #val(days, "freq_bias", sample: "all_products", product: "GPM", stratifier: "all", level: "all", threshold: "0.1") and #val(days, "freq_bias", sample: "all_products", product: "GPM", stratifier: "all", level: "all", threshold: "1") for GPM) — its fewer, heavier false alarms cost it less at the daily scale.

*The longest dry spell* at each gauge is reproduced in the typical case: at 1 mm/day the median over gauges of the longest spell is #dd("median_station_max", "gpm_gsmap_full", "GPM", "Gauge", "1") days at the gauges, #dd("median_station_max", "gpm_gsmap_full", "GPM", "GPM", "1") for GPM and #dd("median_station_max", "gpm_gsmap_full", "GSMaP", "GSMaP", "1") for GSMaP. The gauges' record, #dd("max_length", "gpm_gsmap_full", "GPM", "Gauge", "1") days, against #dd("max_length", "gpm_gsmap_full", "GPM", "GPM", "1") for GPM, is a silent gauge (@sec-silent).

#figure(
  image("figures/norain_fig9_spells.png", width: 100%),
  caption: [Left: the share of dry hours in complete dry spells of at least $L$ hours, full record. Middle: wet-day frequency bias at 0.1 and 1 mm/day. Right: the longest run of days below 1 mm at each gauge, gauge against product.],
) <fig-spells>

== Calibrated and uncalibrated IMERG <sec-calibration>

IMERG Final's gauge calibration is a monthly adjustment of the multi-satellite field toward the GPCC gauge analysis. It does not remove dry-hour rain; it adds to it (@fig-calibration). The uncalibrated field reports ≥ 0.1 mm/h on #ov("POFD", "gpm_cal_uncal", "GPM_uncal", "0.1")% of gauge-dry hours and the calibrated one on #ov("POFD", "gpm_cal_uncal", "GPM", "0.1")% (#ov("POFD", "gpm_cal_uncal", "GPM_uncal", "0.5")% and #ov("POFD", "gpm_cal_uncal", "GPM", "0.5")% at 0.5 mm/h); the spurious rain grows from #ov("spurious_mm", "gpm_cal_uncal", "GPM_uncal", "0.1") to #ov("spurious_mm", "gpm_cal_uncal", "GPM", "0.1") mm per year while the share of the volume on dry hours stays at #ov("dry_volume_share", "gpm_cal_uncal", "GPM", "0.1")%. The calibration scales the amounts (calibrated/uncalibrated volume #val(calibration, "volume_ratio_calibrated_to_uncalibrated", kind: "pattern")) and so lifts some drizzle over the thresholds, but it leaves the occurrence pattern alone: both series are zero on #val(calibration, "both_zero", kind: "pattern")% of cells and disagree about zero on only #val(calibration, "zero_in_calibrated_only", kind: "pattern") + #val(calibration, "zero_in_uncalibrated_only", kind: "pattern")%. In winter the calibration helps: it doubles GPM's wet days at 1 mm (wet-day bias #val(days, "freq_bias", sample: "gpm_cal_uncal", product: "GPM_uncal", stratifier: "season", level: "DJF", threshold: "1") uncalibrated, #val(days, "freq_bias", sample: "gpm_cal_uncal", product: "GPM", stratifier: "season", level: "DJF", threshold: "1") calibrated), where the multi-satellite field misses rain. The two series share their clock: the pooled correlation peaks at lag 0 (r = #val(calibration, "r", kind: "lag", lag: "0"), against #val(calibration, "r", kind: "lag", lag: "-1") and #val(calibration, "r", kind: "lag", lag: "1") at −1 h and +1 h).

The dry-hour rain is therefore made upstream, in the multi-satellite occurrence field, and no monthly gauge calibration can reach it. Only an hourly, local gauge correction such as the fusion methods can.

#figure(
  image("figures/norain_fig10_calibration.png", width: 100%),
  caption: [Calibrated (filled) and uncalibrated (hollow, dashed) IMERG on gauge-dry hours, full record: exceedance, false-alarm rate by season, and the clock check.],
) <fig-calibration>

== Interpolation and fusion on dry hours <sec-fusion>

The benchmark's held-out predictions, scored on the #fu("n_dry", "GPM", "adw") gauge-dry station-hours of the GPM anchor's mask (balanced spatial CV unless stated), are in @tab-fusion and @fig-fusion. Rain per year here is per year of the benchmark's scored hours, which follow FY4B's coverage and are weighted to summer, so it compares methods with each other, not with the full-record rates above.

*Gauge-only interpolation is not dry either.* ADW predicts ≥ 0.1 mm/h on #fu("POFD01", "GPM", "adw")% of the held-out gauges' dry hours — more often than raw GPM (#fu("POFD01", "GPM", "raw")%) — and puts as much rain on them (#fu("spurious_mm", "GPM", "adw") against #fu("spurious_mm", "GPM", "raw") mm per year); it is exactly zero on only #fu("zero_share", "GPM", "adw")% of them, against GPM's #fu("zero_share", "GPM", "raw")%. An interpolator carries the rain of wet training gauges into the dry held-out one, and the farther away they are the more it smears: ADW's dry-hour mean grows from #val(fusion-strata, "mean_dry", anchor: "GPM", method: "adw", stratifier: "distance", level: "0-20 km") mm/h within 20 km of the nearest training gauge to #val(fusion-strata, "mean_dry", anchor: "GPM", method: "adw", stratifier: "distance", level: "50-100 km") at 50–100 km, where raw GPM stays at #val(fusion-strata, "mean_dry", anchor: "GPM", method: "raw", stratifier: "distance", level: "50-100 km"). Under random CV, where a held-out gauge keeps its neighbours, ADW is wet on only #fu("POFD01", "GPM", "adw", scheme: "random")% of dry hours. ADW nevertheless has the lower dry-hour RMSE (#fu("RMSE_dry", "GPM", "adw") against #fu("RMSE_dry", "GPM", "raw") mm/h): its errors are many and small, the satellite's few and large, and squared error punishes the large ones.

*Unblended fusion imports the anchor's false rain.* MGWR on GPM is drier than ADW where the anchor is dry (dry-hour mean #val(fusion-strata, "mean_dry", anchor: "GPM", method: "auto", stratifier: "anchor_state", level: "anchor dry") against #val(fusion-strata, "mean_dry", anchor: "GPM", method: "adw", stratifier: "anchor_state", level: "anchor dry") mm/h) but on the #val(fusion-strata, "n_dry", anchor: "GPM", method: "adw", stratifier: "anchor_state", level: "anchor wet") dry hours where GPM is wet it predicts #val(fusion-strata, "mean_dry", anchor: "GPM", method: "auto", stratifier: "anchor_state", level: "anchor wet") mm/h and is wet #val(fusion-strata, "POFD01", anchor: "GPM", method: "auto", stratifier: "anchor_state", level: "anchor wet")% of the time (ADW #val(fusion-strata, "mean_dry", anchor: "GPM", method: "adw", stratifier: "anchor_state", level: "anchor wet") and #val(fusion-strata, "POFD01", anchor: "GPM", method: "adw", stratifier: "anchor_state", level: "anchor wet")%). The net: its dry-hour RMSE is #fu("dry_vs_adw", "GPM", "mgwr").replace("−", "")% larger than ADW's #ci(fu("dry_vs_adw_ci", "GPM", "mgwr")), and dry hours make up #fu("gap_dry_total", "GPM", "mgwr")% of its whole squared-error gap to ADW (@fig-fusion-where, left) — almost all of it on hours when the gauge is dry and the anchor wet. On FY4B the dry-hour loss is #fu("dry_vs_adw", "FY4B", "mgwr").replace("−", "")%, on GSMaP #fu("dry_vs_adw", "GSMaP", "mgwr").replace("−", "")%. This is the mechanism behind the fusion report's finding that unblended fusion loses to ADW overall.

*The blend repairs it.* Pulling the correction toward ADW where the satellite reports rain cuts GPM's anchor-wet dry-hour mean to #val(fusion-strata, "mean_dry", anchor: "GPM", method: "blend_mgwr", stratifier: "anchor_state", level: "anchor wet") mm/h. The MGWR + ADW blend then beats ADW on the dry hours themselves by #fu("dry_vs_adw", "GPM", "blend_mgwr")% #ci(fu("dry_vs_adw_ci", "GPM", "blend_mgwr")) and puts less rain on them than any other method (#fu("spurious_mm", "GPM", "blend_mgwr") mm per year); the agreement blend gains #fu("dry_vs_adw", "GPM", "blend_agrenv_mgwr")% #ci(fu("dry_vs_adw_ci", "GPM", "blend_agrenv_mgwr")). On FY4B, GSMaP and the equal-weight merge the blends stay a few percent behind ADW on dry hours (#fu("dry_vs_adw", "FY4B", "blend_mgwr")%, #fu("dry_vs_adw", "GSMaP", "blend_mgwr")%, #fu("dry_vs_adw", "MERGED_MEAN", "blend_mgwr")%); whether they beat ADW overall is then decided on the wet hours (see the fusion report). Under random CV the GPM blend ties ADW on dry hours (#fu("dry_vs_adw", "GPM", "blend_mgwr", scheme: "random")% #ci(fu("dry_vs_adw_ci", "GPM", "blend_mgwr", scheme: "random"))).

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.4pt),
    table.header(rule, [*Anchor*], [*Method*], [*Dry mean*], [*≥ 0.1 %*], [*Zero %*], [*Rain mm/yr*], [*Dry RMSE*],
      [*Dry RMSE vs ADW %*], thin),
    ..(("GPM", ("raw", "adw", "idw", "tps", "gwr", "mgwr", "blend_mgwr", "blend_agrenv_mgwr")),
       ("MERGED_OLS", ("raw", "mgwr", "blend_mgwr", "blend_agrenv_mgwr"))).map(((anchor, methods)) => methods.map(m => {
      let r = find(fusion, scheme: "balanced_spatial", anchor: anchor, method: m)
      (A.at(anchor), M.at(m), r.mean_dry, r.POFD01, r.zero_share, r.spurious_mm, r.RMSE_dry,
        if r.dry_vs_adw == "" [—] else [#r.dry_vs_adw.replace("-", "−") #ci(r.dry_vs_adw_ci)])
    })).flatten(),
    rule,
  ),
  caption: [The methods on held-out gauges' dry hours, balanced spatial CV. Dry mean and RMSE in mm/h; ≥ 0.1: share of dry hours predicted ≥ 0.1 mm/h; Zero: share predicted exactly 0; last column: 1 − dry-hour RMSE / ADW's on identical cells (+ = better than ADW), 95% paired day-block bootstrap. Gauge-only rows are the same on every anchor's mask to within rounding.],
) <tab-fusion>
]

*The OLS merge is never dry.* Its fitted intercept leaves a small positive value on almost every hour: the raw merged field is exactly zero on #fu("zero_share", "MERGED_OLS", "raw")% of dry hours, and after MGWR and the blends on #fu("zero_share", "MERGED_OLS", "blend_agrenv_mgwr")–#fu("zero_share", "MERGED_OLS", "mgwr")%. Summed over a day this reaches 0.1 mm everywhere, so every held-out gauge-dry day is predicted wet (daily false-alarm rate #val(fusion-days, "POFD", anchor: "MERGED_OLS", method: "blend_mgwr", stratifier: "all", level: "all", threshold: "0.1")%, wet-day frequency bias #val(fusion-days, "freq_bias", anchor: "MERGED_OLS", method: "blend_mgwr", stratifier: "all", level: "all", threshold: "0.1")). Yet the raw OLS merge has the smallest dry-hour RMSE of any method (#fu("RMSE_dry", "MERGED_OLS", "raw") mm/h, #fu("dry_vs_adw", "MERGED_OLS", "raw")% better than ADW): a thin uniform drizzle is exactly what squared error rewards when the truth is 0 most of the time and uncertain the rest. Part of the OLS merge's RMSE lead in the fusion report is therefore bought with a field that has no dry days.

*Zeroing small values* (@fig-fusion-where, right; descriptive, not tuned in-fold). Setting every prediction below 0.3 mm/h to zero lowers ADW's dry-hour false-alarm rate from #val(fusion-zero, "POFD", anchor: "GPM", method: "adw", tau: "0") to #val(fusion-zero, "POFD", anchor: "GPM", method: "adw", tau: "0.3")% and raises its CSI at 0.1 mm/h from #val(fusion-zero, "CSI", anchor: "GPM", method: "adw", tau: "0") to #val(fusion-zero, "CSI", anchor: "GPM", method: "adw", tau: "0.3"), but makes the overall RMSE slightly worse for every method (ADW #val(fusion-zero, "RMSE_change", anchor: "GPM", method: "adw", tau: "0.3")%, MGWR + ADW blend #val(fusion-zero, "RMSE_change", anchor: "GPM", method: "blend_mgwr", tau: "0.3")%), because the same cut removes light rain on wet hours. A dry gate is a choice for occurrence and dry-day statistics, and an RMSE-tuned model will never make it on its own.

*Silent gauges in the benchmark.* #val(fusion-screened, "n", anchor: "GPM", method: "mgwr", level: "gauge dry, silent-gauge spell") of the held-out dry station-hours fall inside flagged silent-gauge spells. On them every method predicts about twice its usual dry-hour amount (ADW #val(fusion-screened, "mean_ref", anchor: "GPM", method: "mgwr", level: "gauge dry, silent-gauge spell") against #val(fusion-screened, "mean_ref", anchor: "GPM", method: "mgwr", level: "gauge dry") mm/h) — the neighbours were wet — and is scored as wrong for it. They account for only #val(fusion-screened, "gap_share", anchor: "GPM", method: "mgwr", level: "gauge dry, silent-gauge spell")% of MGWR's squared-error gap to ADW, so they do not change any ranking.

#figure(
  image("figures/norain_fig11_fusion_dry.png", width: 100%),
  caption: [What each method predicts on held-out gauges' dry hours, by anchor, balanced spatial CV; gauge-only methods in ink.],
) <fig-fusion>

#figure(
  image("figures/norain_fig12_fusion_where.png", width: 100%),
  caption: [GPM anchor, balanced spatial CV. Left: each gauge × anchor quadrant's contribution to the gap between the method's squared error and ADW's. Middle: dry-hour mean by distance to the nearest training gauge and to the gauge's own rain. Right: overall RMSE change when predictions below τ are set to zero.],
) <fig-fusion-where>

= Answers

+ *How often, how much, when?* On #ov("POFD", "gpm_gsmap_full", "GPM", "0.1")–#ov("POFD", "gpm_gsmap_full", "GSMaP", "0.1")% of gauge-dry hours at 0.1 mm/h (GPM, GSMaP; FY4B #ov("POFD", "all_products", "FY4B", "0.1")% on its hours), #ov("POFD", "gpm_gsmap_full", "GPM", "0.5")–#ov("POFD", "gpm_gsmap_full", "GSMaP", "0.5")% at 0.5 mm/h, mostly drizzle; a third of GPM's and GSMaP's volume and #ov("dry_volume_share", "all_products", "FY4B", "0.1")% of FY4B's. Most in summer and in the afternoon (GPM) or late morning (GSMaP, spuriously), uniformly across gauges and relief.
+ *Is it really false?* Mostly not in the sense of "rain from nothing". The rate rises twenty-fold next to the gauge's own rain and next to wet neighbours, falls to #st("POFD05", "gpm_gsmap_full", "GPM", "network", "network dry")% on hours when the whole network is dry, and at 0.5 mm/h sits inside the range of a neighbouring gauge. What remains is drizzle below the bucket's resolution, rain evaporating below the cloud, and — for FY4B — genuine phantom rain. Silent gauges explain a few percent of the false alarms.
+ *Dry spells and dry days?* Drizzle breaks the long hourly dry spells at 0.1 mm/h, and the products have too few dry days (#dd("dry_share", "gpm_gsmap_full", "GPM", "GPM", "0.1")% for GPM against #dd("dry_share", "gpm_gsmap_full", "GPM", "Gauge", "0.1")% at the gauges). At 0.5 mm/h and 1 mm/day the dry-spell structure and the typical longest dry spell are reproduced. IMERG's calibration does not change any of it.
+ *The fusion methods?* All of them rain on held-out dry hours, gauge-only interpolation as often as the satellite. Unblended residual correction inherits the anchor's false rain, which is why it loses to ADW; blending toward ADW where the satellite is wet removes that loss. RMSE rewards thin drizzle over sharp zeros, which the OLS merge exploits.

= What this means for the correction

- *Blend, and keep blending.* The dry hours are where unblended fusion loses; the ADW blend is the part of the method that fixes them. Any new anchor or correction model should be judged on its anchor-wet dry hours first.
- *Treat the satellite's wet signal as a probability of nearby rain, not of rain at the gauge.* A product wet next to a gauge's rain, or wet where the others are too, is informative; the same value on a network-dry hour, or in dry air, is not. Product agreement and ERA5-Land dew-point depression are both available at prediction time and could enter the blend weight or a dry gate.
- *Do not expect a monthly calibration, or an average of products, to remove dry-hour rain.* Both leave the occurrence pattern as it is.
- *Look at more than RMSE.* The OLS merge wins dry-hour RMSE with a field that has no dry days. If dry-day frequency or dry spells matter for the application (drought indices, runoff), the benchmark needs an occurrence score — dry-day frequency bias, POFD at 0.1 mm/h — beside RMSE when it selects methods.
- *Screen the gauges.* Silent-gauge spells should be masked before training and scoring; the rule of @sec-silent-method is cheap and its flags are in `silent_gauge_spells.csv`.

= Limitations

- *Point against pixel.* Every "false alarm" is a comparison of a 200 cm² funnel with a 10 km (GPM, GSMaP) or 4 km (FY4B) pixel. The gauge-to-gauge baseline bounds this but cannot remove it.
- *The gauges' resolution.* With a 0.5 mm bucket, an hour with 0.1–0.4 mm of rain can read 0.0; some of the drizzle counted as false is real. The 0.5 mm/h results are the fairer ones.
- *Spurious mm per year* is only meaningful on the full record. The FY4B-hours sample and the benchmark grid are weighted to summer, so there the share of the volume on dry hours is used instead, and rain per year compares methods only.
- *Partial days and censored spells.* Hourly dry spells exist only on the full record: FY4B never has 23:00, so every spell is cut at under a day. Daily results on the FY4B hours and the benchmark grid use days with at least 20 of 24 hours scored, summed over the same hours for gauge and product.
- *Winter.* Tipping buckets do not record snow as it falls, so some winter "false alarms" are real snowfall; the silent-gauge screen catches the longest of these spells but not single days.
- *The silent-gauge screen is a description,* tuned on one rate (#val(silent-sensitivity, "p_dry_given_wet", min_wet_days: "3")% of neighbour-wet days dry at a working gauge) and not validated against station logs. No gauge value was changed.
- *The zero-threshold sensitivity is descriptive.* The threshold was not tuned inside the folds, so it shows what a dry gate could do, not what it would do.

= Next steps

+ *An in-fold dry gate.* Tune a zero threshold, or run the disabled `hurdle_gwr` model, inside the benchmark's folds and score it on RMSE and on occurrence together. This needs a new benchmark run: about 3.5 h for spatial CV only, 7 h for both schemes.
+ *Mask silent gauges in the benchmark,* in training and in scoring, and re-run the claim assessment to see whether any ranking moves.
+ *Condition the blend on what predicts false rain:* product agreement, distance to the nearest wet gauge, ERA5-Land dew-point depression.
+ *Add an occurrence criterion to method selection* (dry-day frequency bias, POFD at 0.1 mm/h), and report the OLS merge's dry-day behaviour next to its RMSE in the fusion report.
+ *Check the silent-gauge spells against station maintenance logs,* if they can be obtained.

#pagebreak()

#heading(numbering: none)[Appendix — Reproduction]

From the repository root, one Julia job at a time (the two drivers each hold about 1.3 GB):

```sh
julia -t 4 --project=. scripts/run_no_rain_evaluation.jl         # ~15 min
julia -t 4 --project=. scripts/run_no_rain_fusion_evaluation.jl  # ~35 min
py -3.13 scripts/plot_no_rain_evaluation.py
typst compile Interim_results/no_rain_evaluation.typ
```

The drivers write to `output/no_rain_evaluation/` and its `fusion/` subfolder; the plot script writes the figures there and copies them, with the report tables `tables/norain_*.json`, into this folder. The analysis is in `src/benchmark/NoRainEvaluation.jl` (tests: `test/test-no-rain-evaluation.jl`). Two gates run before anything is written: the satellite driver checks every gauge-dry count and rate against the temporal evaluation's `false_alarm_metrics.csv`, and the fusion driver checks every method's dry-hour n, RMSE and bias against the benchmark run's own `metrics_pooled.csv`.

- Benchmark run read by the fusion part: #text(size: 8.5pt, raw-val(provenance, "value", key: "run_directory").replace("_", "_\u{200B}")), written at commit #raw(raw-val(provenance, "value", key: "git_commit").slice(0, 7)) on #raw(raw-val(provenance, "value", key: "git_branch")) (working tree dirty: #raw-val(provenance, "value", key: "git_dirty")).
- Bootstrap: #raw-val(provenance, "value", key: "setting_bootstrap_reps") replicates, blocks of one #raw-val(provenance, "value", key: "setting_bootstrap_block"), seed #raw-val(provenance, "value", key: "setting_seed").
- Silent-gauge rule: #raw(raw-val(provenance, "value", key: "setting_silent_gauge_rule")).
- ERA5-Land and uncalibrated IMERG placed on the gauge clock with a +#raw-val(provenance, "value", key: "setting_era5_offset_hours") h offset (UTC to BJT hour-ending).
