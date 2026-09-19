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
  auto: "GWR family, chosen in-fold", residual_gwr: "Residual GWR", mixed_gwr: "Mixed GWR", mgwr: "MGWR",
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

KEY-FINDINGS

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

RECORD

== How often the products rain on dry hours <sec-how-often>

HOW-OFTEN

== When and where <sec-when-where>

WHEN-WHERE

== Is the rain really false? <sec-false>

FALSE

== Silent gauges <sec-silent>

SILENT

== Dry spells and dry days <sec-spells>

SPELLS

== Calibrated and uncalibrated IMERG <sec-calibration>

CALIBRATION

== Interpolation and fusion on dry hours <sec-fusion>

FUSION

= Answers

ANSWERS

= What this means for the correction

IMPLICATIONS

= Limitations

LIMITATIONS

= Next steps

NEXT

#pagebreak()

#heading(numbering: none)[Appendix — Reproduction]

APPENDIX
