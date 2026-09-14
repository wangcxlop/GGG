// Interim results: what combining rain gauges with FY4B, GPM and GSMaP gains, which interpolation and
// fusion methods do best, and which product - or merge of products - to build on.
// Compile from this folder:  typst compile gauge_satellite_fusion.typ
// Every number is read from tables/*.json, which scripts/plot_gauge_satellite_fusion.py writes from
// output/gauge_satellite_fusion/*.csv; figures are copies of output/gauge_satellite_fusion/figures/.

#set document(
  title: "Gauges and satellites: interpolation, fusion and product choice",
  date: datetime(year: 2026, month: 9, day: 14),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Gauge–satellite fusion — interim results
      #h(1fr) 14 September 2026
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
// Data. Tables are loaded once; the prose reads single values through `val`, so a number in the text
// is always the number in the table.

#let load(name) = json("tables/" + name + ".json")
#let overall = load("overall_metrics")
#let ranking = load("interpolation_ranking")
#let fusion = load("fusion_vs_adw")
#let anchors = load("anchor_metrics")
#let versus-gpm = load("anchor_vs_gpm")
#let detection = load("detection_skill")
#let weights = load("merge_weights")
#let blends = load("blend_weights")
#let provenance = load("provenance")
#let robust = load("robustness")
#let where = load("where_it_helps")

#let find(rows, ..conditions) = {
  let hits = rows.filter(row => conditions.named().pairs().all(pair => row.at(pair.at(0)) == pair.at(1)))
  if hits.len() == 0 { none } else { hits.at(0) }
}
/// One cell as typeset text: ASCII hyphens in numbers become minus signs.
#let val(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else { hit.at(column).replace("-", "−") }
}
/// A percentage for prose. A leading plus is dropped - the sentence already says "beats" or "falls by" -
/// while a minus sign is kept.
#let pct(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else [#hit.at(column).trim("+", at: start).replace("-", "−")%]
}
/// The smallest and largest value of `column` over the rows matching `keep`, as "a–b" (signs dropped).
#let span(rows, column, keep) = {
  let values = rows.filter(keep).map(row => float(row.at(column).replace("+", "")))
  if values.len() == 0 { text(fill: red)[?] } else {
    let low = calc.min(..values)
    let high = calc.max(..values)
    [#calc.round(low, digits: 1)–#calc.round(high, digits: 1)]
  }
}

#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)
#let ci(body) = text(size: 6.5pt, fill: luma(90), body.replace("-", "−"))
#let note(body) = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, fill: luma(246),
  stroke: (left: 2.5pt + rgb("#2a78d6")), body)

#let M = (
  raw: "Raw satellite", adw: "ADW (gauges only)", idw: "IDW", tps: "TPS", gwr: "GWR (gauges only)",
  infold: "GWR family, chosen in-fold", rgwr: "Residual GWR", mgwr: "MGWR", bmgwr: "MGWR + ADW blend",
  amgwr: "MGWR + agreement blend",
)
#let anchor-label = (FY4B: "FY4B", GPM: "GPM", GSMaP: "GSMaP", MERGED_MEAN: "Merged (equal)",
  MERGED_OLS: "Merged (OLS)", "Gauges only": "Gauges only")

// ---------------------------------------------------------------------------------------------
// Title

#align(left)[
  #set text(hyphenate: false)
  #text(font: ("Segoe UI", "Arial"), size: 9pt, fill: luma(100), tracking: 0.04em)[INTERIM RESULTS · 14 SEPTEMBER 2026]
  #v(0.2em)
  #text(font: ("Segoe UI", "Arial"), size: 19pt, weight: "semibold")[What do rain gauges and satellites gain from each other?]
  #v(0.1em)
  #text(size: 12pt, fill: luma(60))[Interpolation versus gauge–satellite fusion at 237 hourly gauges, 2022–2024, with FY4B, GPM IMERG, GSMaP and two merges of the three]
]
#v(0.6em)
#line(length: 100%, stroke: 0.5pt + luma(180))

= Key findings

- *Gauges repair every satellite product.* At gauges held out of the fit, correcting a product with the surrounding gauges cuts its hourly RMSE by #pct(overall, "vs_raw", anchor: "FY4B", method: M.amgwr) for FY4B, #pct(overall, "vs_raw", anchor: "GPM", method: M.amgwr) for GPM and #pct(overall, "vs_raw", anchor: "GSMaP", method: M.amgwr) for GSMaP (MGWR + agreement blend), and raises the correlation from #val(overall, "r", anchor: "FY4B", method: M.raw), #val(overall, "r", anchor: "GPM", method: M.raw) and #val(overall, "r", anchor: "GSMaP", method: M.raw) to #val(overall, "r", anchor: "FY4B", method: M.amgwr), #val(overall, "r", anchor: "GPM", method: M.amgwr) and #val(overall, "r", anchor: "GSMaP", method: M.amgwr).
- *Satellites add to the gauges only in the right fusion.* Gauge-only ADW (RMSE #val(overall, "RMSE", anchor: "GPM", method: M.adw)) beats every plain residual-correction model on a single product. Blending the correction back toward ADW where the satellite reports rain turns that around: #M.amgwr beats ADW by #pct(fusion, "change", stratum: "overall", anchor: "GPM", method: M.amgwr) on GPM #ci(val(fusion, "ci", stratum: "overall", anchor: "GPM", method: M.amgwr)) and by #pct(fusion, "change", stratum: "overall", anchor: "MERGED_OLS", method: M.amgwr) on the OLS merge #ci(val(fusion, "ci", stratum: "overall", anchor: "MERGED_OLS", method: M.amgwr)).
- *The satellite pays where the gauges cannot see.* The blended MGWR's advantage over ADW grows with distance from the nearest training gauge — on GPM from #pct(where, "change", anchor: "GPM", method: M.amgwr, group: "nearest_train_km", level: "0_20") within 20 km to #pct(where, "change", anchor: "GPM", method: M.amgwr, group: "nearest_train_km", level: "50_100") at 50–100 km. Unblended fusion is best at detecting ≥ 8 mm/h rain (CSI #val(detection, "CSI", anchor: "GPM", method: M.infold, threshold: "8") against ADW's #val(detection, "CSI_adw", anchor: "GPM", method: M.infold, threshold: "8") on GPM) but loses badly on dry hours, which is why it loses overall.
- *Best interpolator: ADW*, statistically ahead of IDW by a small margin (#pct(ranking, "change", stratum: "overall", method: M.idw)) and clearly ahead of TPS (#pct(ranking, "change", stratum: "overall", method: M.tps)) and coordinate GWR (#pct(ranking, "change", stratum: "overall", method: M.gwr)).
- *Best fusion method: MGWR blended toward ADW.* Of the three GWR variants MGWR is the strongest on every anchor, and the two blending rules perform about the same; whether to blend matters more than which variant or rule.
- *Best single product: GPM; better still, merge all three.* On identical cells GPM's raw field (RMSE #val(anchors, "RMSE", anchor: "GPM", method: M.raw)) beats GSMaP (#val(anchors, "RMSE", anchor: "GSMaP", method: M.raw)) and FY4B (#val(anchors, "RMSE", anchor: "FY4B", method: M.raw)), and an OLS merge fitted on training gauges beats GPM by #pct(versus-gpm, "change", anchor: "MERGED_OLS", method: M.raw) as a raw field and by #pct(versus-gpm, "change", anchor: "MERGED_OLS", method: M.mgwr) after MGWR. Once blended the merge's lead shrinks to #pct(versus-gpm, "change", anchor: "MERGED_OLS", method: M.amgwr).

= The questions

The supervisor asked three things, and they structure this report:

+ What does integrating gauge observations with satellite precipitation gain? Show it with metrics and figures (@sec-benefit, @sec-where).
+ Which interpolation method performs best, and which data-fusion method (@sec-interpolation, @sec-fusion)?
+ Which satellite product has the highest quality — or should all three be integrated (@sec-products)?

= Data and design

- *Gauges.* 237 hourly tipping-bucket gauges in the Shiyan study area (0.5 mm resolution).
- *Products.* FY4B precipitation (strict hourly aggregation, navigation-corrected), GPM IMERG V07 and GSMaP V8 operational, each sampled at the gauge pixel and aligned to the gauge clock.
- *Period.* January 2022 – December 2024, restricted to the *13,471 hours* on which all three products are available. FY4B's strict completeness sets this grid, so the comparison is between products on identical hours.
- *Validation.* Five-fold cross-validation by *station*. The primary scheme holds out spatially compact groups of gauges (balanced spatial folds), so a held-out gauge is typically ~24 km from the nearest training gauge — the situation of an ungauged area. A random station split, where a held-out gauge keeps its neighbours, is reported as a secondary check (@sec-robustness).
- *No selection sees a held-out gauge.* Every hyperparameter (kernel, bandwidth, IDW power and neighbours, TPS smoothing, the residual shrinkage and the blending weights) is tuned on an inner spatial split of the training gauges. Covariates (ERA5-Land, MODIS NDVI and DEM terrain) are selected inside every training fold by permutation tests. The OLS merge weights are also fitted inside every fold.
- *Metrics.* RMSE, MAE, bias and Pearson r over all station-hours; POD, FAR and CSI at 0.1, 2.5, 8 and 16 mm/h (hit when both gauge and estimate reach the threshold). Rain-intensity strata follow the gauge: dry < 0.1, light 0.1–2.5, moderate 2.5–8, heavy ≥ 8 mm/h. These are the benchmark's strata and differ from the temporal-evaluation report's intensity classes.
- *Uncertainty.* Paired day-block bootstrap (2,000 replicates): whole days are resampled for both methods together. P-values are Holm-corrected across the whole family of comparisons in a table, not only within one row.

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    stroke: none,
    inset: (x: 5pt, y: 2.6pt),
    table.header(rule, [*Family*], [*Method*], [*What it does*], thin),
    [Satellite only], [Raw], [the product's value at the gauge pixel],
    table.cell(rowspan: 4)[Gauges only], [IDW], [inverse-distance weighting of the training gauges],
    [ADW], [angular distance weighting: IDW with a penalty for gauges that lie in the same direction],
    [TPS], [thin-plate spline],
    [GWR], [geographically weighted regression on the coordinates alone],
    table.cell(rowspan: 4)[Fusion], [Residual GWR], [models the gauge − satellite residual with GWR on coordinates and the fold's covariates, then adds it back to the satellite field],
    [Mixed GWR], [the same, with some covariates held global],
    [MGWR], [multiscale GWR: each term gets its own bandwidth],
    [GWR family, chosen in-fold], [one of the four GWR models, picked per fold on the inner split — the selection-free answer],
    table.cell(rowspan: 2)[Blended fusion], [\+ ADW blend], [where the satellite reports ≥ 0.1 mm/h, $(1 - lambda) dot "fusion" + lambda dot "ADW"$, $lambda$ tuned per fold],
    [\+ agreement blend], [one $lambda$ per band of (how many of the three products report rain) × (their maximum)],
    table.cell(rowspan: 2)[Merged anchors], [Merged (equal)], [mean of FY4B, GPM and GSMaP],
    [Merged (OLS)], [gauge ≈ $a + b_1 "FY4B" + b_2 "GPM" + b_3 "GSMaP"$, fitted on the fold's training gauges],
    rule,
  ),
  caption: [Methods. Every fusion and blend method runs on each of the five anchors.],
) <tab-methods>
]

#note[
  *Two things to keep in mind when reading "best".* Picking the top cell of a table of about fifty configurations is a selection made on held-out scores. The in-fold choice ("GWR family, chosen in-fold") is reported beside it as the answer that involves no such selection, and the Holm correction covers the whole family. Second, the agreement blend reads all three products to decide where to blend, so "MGWR + agreement blend on GPM" already uses FY4B and GSMaP.
]

= Results

== What integration gains <sec-benefit>

@fig-benefit and @tab-benefit put four stages side by side for each anchor: the raw satellite field, gauges alone (ADW), fusion chosen inside the fold, and blended MGWR. Adding gauges to a satellite product helps every product, and the weaker the product the more it helps: FY4B's RMSE falls by #pct(overall, "vs_raw", anchor: "FY4B", method: M.amgwr), GSMaP's by #pct(overall, "vs_raw", anchor: "GSMaP", method: M.amgwr) and GPM's by #pct(overall, "vs_raw", anchor: "GPM", method: M.amgwr). Bias moves to near zero (@tab-benefit).

Whether the satellite adds to the gauges is a different question and the answer depends on the method. Plain residual correction on a single product loses to ADW: MGWR on GPM is #pct(fusion, "change", stratum: "overall", anchor: "GPM", method: M.mgwr) against ADW, and on FY4B #pct(fusion, "change", stratum: "overall", anchor: "FY4B", method: M.mgwr). The blend reverses this on GPM, GSMaP and both merges.

#figure(
  image("figures/fusion_fig1_benefit.png", width: 100%),
  caption: [Raw satellite, gauges only, and two levels of fusion, for each anchor. Balanced spatial CV.],
) <fig-benefit>

#let benefit-methods = (M.raw, M.adw, M.infold, M.mgwr, M.bmgwr, M.amgwr)
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.4pt),
    table.header(rule, [*Anchor*], [*Method*], [*RMSE*], [*MAE*], [*Bias*], [*r*], [*CSI 2.5*], [*CSI 8*],
      [*vs raw*], [*vs ADW*], thin),
    ..for anchor in ("FY4B", "GPM", "GSMaP", "MERGED_MEAN", "MERGED_OLS") {
      let rows = overall.filter(row => row.anchor == anchor and row.method in benefit-methods)
      (
        table.cell(rowspan: rows.len(), anchor-label.at(anchor)),
        ..for row in rows {
          (row.method, row.RMSE, row.MAE, row.Bias.replace("-", "−"), row.r, row.at("CSI2.5"), row.CSI8,
            row.vs_raw.replace("-", "−") + "%", row.vs_adw.replace("-", "−") + "%")
        },
        thin,
      )
    },
  ),
  caption: [Overall skill at held-out gauges (balanced spatial CV; each anchor on its own evaluation mask). "vs raw" and "vs ADW": relative RMSE improvement.],
) <tab-benefit>

== Best interpolation method <sec-interpolation>

Among the gauge-only methods ADW is best overall (@fig-interpolation, @tab-interpolation). IDW is almost the same estimator and trails it by #pct(ranking, "change", stratum: "overall", method: M.idw), a difference that is significant only because 3 million paired cells make it so. TPS (#pct(ranking, "change", stratum: "overall", method: M.tps)) and GWR on coordinates (#pct(ranking, "change", stratum: "overall", method: M.gwr)) are clearly worse. TPS is the one gauge-only method that does better than ADW on heavy hours (#pct(ranking, "change", stratum: "heavy", method: M.tps)): a spline smooths less than a weighted average and keeps more of a storm's peak.

#figure(
  image("figures/fusion_fig2_interpolators.png", width: 100%),
  caption: [Gauge-only interpolators against ADW, and their detection skill.],
) <fig-interpolation>

#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.4pt),
    table.header(rule, [*Hours*], [*Method*], [*RMSE*], [*ADW*], [*Change*], [*95% CI*], [*Holm p*], thin),
    ..ranking.map(row => (row.stratum, row.method, row.RMSE, row.RMSE_adw, row.change.replace("-", "−") + "%",
      ci(row.ci), row.p_holm)).flatten(),
    rule,
  ),
  caption: [RMSE of each interpolator against ADW on the same cells. Change: positive means better than ADW.],
) <tab-interpolation>
]

== Best fusion method <sec-fusion>

@fig-fusion shows every fusion method on every anchor against ADW. Three patterns hold throughout:

- *MGWR is the strongest GWR variant* on every anchor. Residual and mixed GWR are nearly identical: in 32 of the 50 fold cells the covariate screening assigned no global role, which makes the two the same model (`run_status.csv`, column `duplicate_of`).
- *Unblended fusion wins on heavy hours but loses overall* on the single products. On heavy hours every configuration beats ADW, by up to #pct(fusion, "change", stratum: "heavy", anchor: "MERGED_MEAN", method: M.mgwr) (MGWR on the equal merge). Over all hours the same models lose, and @sec-where shows where: on dry and light hours, where a satellite false alarm enters the correction as a fixed offset.
- *The blend is what makes fusion pay overall.* Leaning toward ADW where the satellite reports rain removes most of the dry-hour loss, at the price of part of the heavy-hour gain: on GPM, MGWR's heavy-hour improvement of #pct(fusion, "change", stratum: "heavy", anchor: "GPM", method: M.mgwr) becomes #pct(fusion, "change", stratum: "heavy", anchor: "GPM", method: M.bmgwr) once blended.

@tab-fusion-top lists the strongest configurations. Their confidence intervals overlap widely, so the order among them is not meaningful; what is clear is that all of them are ahead of ADW.

#figure(
  image("figures/fusion_fig3_methods_by_anchor.png", width: 100%),
  caption: [RMSE improvement over ADW for every fusion method and anchor, all hours and heavy hours.],
) <fig-fusion>

#let top-overall = fusion.filter(row => row.stratum == "overall").slice(0, 10)
#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.4pt),
    table.header(rule, [*Anchor*], [*Method*], [*RMSE*], [*ADW*], [*Change*], [*95% CI*], [*Holm p*], thin),
    ..top-overall.map(row => (anchor-label.at(row.anchor), row.method, row.RMSE, row.RMSE_adw,
      row.change.replace("-", "−") + "%", ci(row.ci), row.p_holm_family)).flatten(),
    rule,
  ),
  caption: [The ten fusion configurations with the largest overall RMSE improvement over ADW. Holm p across all 45 configurations.],
) <tab-fusion-top>
]

The blending weights are stable across folds (@tab-blends): on the single products the constant blend moves about two-thirds of the way toward ADW on satellite-wet cells, and on the merges well under half — the merged anchor needs less repair.

#block(breakable: false)[
#figure(
  table(
    columns: (1fr, auto, auto, auto),
    align: (left, right, right, right),
    stroke: none,
    inset: (x: 6pt, y: 2.4pt),
    table.header(rule, [*Anchor*], [*Mean λ*], [*Range over folds*], [*Inner-split RMSE gain*], thin),
    ..blends.map(row => (anchor-label.at(row.anchor), row.lambda_mean, row.lambda_range,
      row.inner_gain.replace("-", "−") + "%")).flatten(),
    rule,
  ),
  caption: [Weight of ADW in the constant blend (MGWR + ADW blend), chosen per fold on the inner split.],
) <tab-blends>
]

== Where the satellite adds to the gauges <sec-where>

@fig-where breaks the comparison with ADW down.

- *By intensity.* The in-fold GWR choice loses heavily on dry hours (GPM: #pct(where, "change", anchor: "GPM", method: M.infold, group: "rain_intensity", level: "no_rain")) and light hours (#pct(where, "change", anchor: "GPM", method: M.infold, group: "rain_intensity", level: "light")) and gains on heavy hours (#pct(where, "change", anchor: "GPM", method: M.infold, group: "rain_intensity", level: "heavy")). The blend removes the dry-hour loss (#pct(where, "change", anchor: "GPM", method: M.amgwr, group: "rain_intensity", level: "no_rain")) and keeps a smaller heavy-hour gain (#pct(where, "change", anchor: "GPM", method: M.amgwr, group: "rain_intensity", level: "heavy")).
- *By distance.* The blend's advantage over ADW grows with the distance from the held-out gauge to the nearest training gauge on every anchor; on the OLS merge from #pct(where, "change", anchor: "MERGED_OLS", method: M.amgwr, group: "nearest_train_km", level: "0_20") within 20 km to #pct(where, "change", anchor: "MERGED_OLS", method: M.amgwr, group: "nearest_train_km", level: "50_100") at 50–100 km. Close to other gauges interpolation is hard to beat; far from them the satellite field carries information the network does not.
- *Detection.* At ≥ 8 mm/h raw GPM and GSMaP detect as well as ADW (ΔCSI #val(detection, "delta", anchor: "GPM", method: M.raw, threshold: "8") and #val(detection, "delta", anchor: "GSMaP", method: M.raw, threshold: "8")), raw FY4B worse (#val(detection, "delta", anchor: "FY4B", method: M.raw, threshold: "8")). The in-fold GWR choice improves on ADW most clearly there (@tab-detection; ΔCSI #val(detection, "delta", anchor: "GPM", method: M.infold, threshold: "8") on GPM, #val(detection, "delta", anchor: "GSMaP", method: M.infold, threshold: "8") on GSMaP); the blend keeps only part of that gain (#val(detection, "delta", anchor: "GPM", method: M.amgwr, threshold: "8") and #val(detection, "delta", anchor: "GSMaP", method: M.amgwr, threshold: "8")) in exchange for its overall RMSE gain.

#figure(
  image("figures/fusion_fig4_where_it_helps.png", width: 100%),
  caption: [RMSE change against ADW by rain intensity and by distance to the nearest training gauge, and CSI change against ADW by threshold.],
) <fig-where>

#let detection-rows = detection.filter(row => row.threshold in ("8",) and row.anchor in ("FY4B", "GPM", "GSMaP", "MERGED_OLS"))
#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.4pt),
    table.header(rule, [*Anchor*], [*Method*], [*POD*], [*FAR*], [*CSI*], [*CSI ADW*], [*ΔCSI*], [*95% CI*], thin),
    ..detection-rows.map(row => (anchor-label.at(row.anchor), row.method, row.POD, row.FAR, row.CSI, row.CSI_adw,
      row.delta.replace("-", "−"), ci(row.ci))).flatten(),
    rule,
  ),
  caption: [Detection of hours ≥ 8 mm/h. ΔCSI: method minus ADW, with a day-block bootstrap interval.],
) <tab-detection>
]

== Which product, or all three <sec-products>

The products are compared on the cells where every anchor was scored (@fig-products, @tab-anchors), with the same method on both sides of each comparison (@tab-versus-gpm).

- *Raw quality: GPM > GSMaP > FY4B.* GPM's field has the lowest RMSE and the highest correlation; FY4B is far behind. This agrees with the temporal evaluation report.
- *After fusion the ranking of the single products holds, but the gaps close.* The blend brings FY4B within #pct(versus-gpm, "change", anchor: "FY4B", method: M.amgwr) of GPM, against #pct(versus-gpm, "change", anchor: "FY4B", method: M.raw) for the raw fields.
- *Merging all three is the strongest anchor.* The OLS merge beats GPM by #pct(versus-gpm, "change", anchor: "MERGED_OLS", method: M.raw) as a raw field and by #pct(versus-gpm, "change", anchor: "MERGED_OLS", method: M.mgwr) after MGWR. Its weights (@tab-weights) put most of the weight on GPM, less on GSMaP and almost none on FY4B, and they barely move between folds. The equal-weight merge is weaker as a raw field but better for heavy hours.

#figure(
  image("figures/fusion_fig5_products.png", width: 100%),
  caption: [Anchors on identical cells, each anchor against GPM, and the OLS merge weights.],
) <fig-products>

#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto),
    align: (left, left, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.4pt),
    table.header(rule, [*Anchor*], [*Method*], [*RMSE*], [*Bias*], [*r*], [*CSI 8*], thin),
    ..anchors.map(row => (anchor-label.at(row.anchor), row.method, row.RMSE, row.Bias.replace("-", "−"), row.r,
      row.CSI8)).flatten(),
    rule,
  ),
  caption: [Every anchor on the cells all five share. "Gauges only": ADW on the same cells.],
) <tab-anchors>
]

#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.4pt),
    table.header(rule, [*Anchor*], [*Method*], [*RMSE*], [*GPM*], [*Change*], [*95% CI*], [*Holm p*], thin),
    ..versus-gpm.map(row => (anchor-label.at(row.anchor), row.method, row.RMSE, row.RMSE_gpm,
      row.change.replace("-", "−") + "%", ci(row.ci), row.p_holm)).flatten(),
    rule,
  ),
  caption: [Each anchor against GPM with the same method, on identical cells. Change: positive means better than GPM.],
) <tab-versus-gpm>
]

#block(breakable: false)[
#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 5pt, y: 2.4pt),
    table.header(rule, [*CV scheme*], [*Intercept*], [*FY4B*], [*GPM*], [*GSMaP*], [*GPM range*], [*Folds fell back*], thin),
    ..weights.map(row => (row.scheme.replace("_", " "), row.intercept, row.FY4B, row.GPM, row.GSMaP, row.GPM_range,
      row.fell_back)).flatten(),
    rule,
  ),
  caption: [OLS merge weights, mean over the five training folds (intercept in mm/h).],
) <tab-weights>
]

== Robustness <sec-robustness>

Under random CV (@fig-robustness, @tab-robustness) a held-out gauge keeps its neighbours and interpolation becomes much easier: ADW's RMSE falls from #val(robust, "RMSE_adw", scheme: "balanced spatial", anchor: "GPM") to #val(robust, "RMSE_adw", scheme: "random", anchor: "GPM"), the in-fold GWR choice is clearly worse than ADW on GPM (RMSE change #pct(robust, "infold_vs_adw", scheme: "random", anchor: "GPM")), and the blend only draws level with it (#pct(robust, "blend_vs_adw", scheme: "random", anchor: "GPM")). The spatial scheme is the relevant estimate for areas without gauges, and there the blend's advantage is present in each year (@tab-robustness, right-hand columns).

#figure(
  image("figures/fusion_fig6_robustness.png", width: 100%),
  caption: [RMSE under spatial and random CV, and the blend's improvement over ADW year by year.],
) <fig-robustness>

#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.4pt),
    table.header(rule, table.cell(rowspan: 2, align: bottom)[*CV*], table.cell(rowspan: 2, align: bottom)[*Anchor*],
      table.cell(rowspan: 2, align: bottom)[*ADW*], table.cell(colspan: 2)[*GWR family, in-fold*],
      table.cell(colspan: 2)[*MGWR + agreement blend*], table.cell(colspan: 3)[*Blend vs ADW by year*],
      [*RMSE*], [*vs ADW*], [*RMSE*], [*vs ADW*], [*2022*], [*2023*], [*2024*], thin),
    ..robust.map(row => (row.scheme, anchor-label.at(row.anchor), row.RMSE_adw, row.RMSE_infold,
      row.infold_vs_adw.replace("-", "−") + "%", row.RMSE_blend, row.blend_vs_adw.replace("-", "−") + "%",
      row.blend_2022.replace("-", "−"), row.blend_2023.replace("-", "−"), row.blend_2024.replace("-", "−"))).flatten(),
    rule,
  ),
  caption: [RMSE under both cross-validation schemes, and the blend's RMSE improvement over ADW (%) by year under spatial CV.],
) <tab-robustness>
]

= Answers

+ *Benefit of integration.* Against the satellite alone the gain is large and certain: #span(overall, "vs_raw", row => row.method == M.amgwr and row.anchor in ("FY4B", "GPM", "GSMaP"))% lower RMSE for the three products with the blended MGWR, much higher correlation and near-zero bias. Against gauges alone it is smaller but real where it matters: the blended MGWR beats ADW by #span(fusion, "change", row => row.stratum == "overall" and row.method == M.amgwr and row.anchor != "FY4B")% overall on GPM, GSMaP and both merges, and by more the further a location is from the gauge network. Unblended fusion gains most on heavy hours (up to #pct(fusion, "change", stratum: "heavy", anchor: "MERGED_MEAN", method: M.mgwr)) but loses over all hours. Close to dense gauges, interpolation alone is about as good.
+ *Best methods.* Interpolation: ADW. Fusion: MGWR residual correction blended toward ADW where the satellite reports rain. Among unblended models MGWR is best, but without the blend no single-product fusion beats ADW overall.
+ *Best product.* GPM is the best single product and FY4B the weakest. Integrating all three through a gauge-fitted OLS merge gives the strongest anchor; the advantage after blending is small but significant. If only one product can be used, use GPM.

= Limitations

- *Hourly, point scale.* Every score is at gauges; 93% of station-hours are dry, so RMSE mostly measures dry-hour behaviour. Daily or areal scores may rank the methods differently.
- *Selection among many configurations.* The Holm correction limits false positives, but the top of @tab-fusion-top is still a selected maximum.
- *The agreement blend already reads all three products*, so on a single-product anchor it is not a single-product method.
- *Covariates.* An earlier exploratory run without ERA5/NDVI/DEM covariates scored the GWR family 0.6–2.2% better. It predates the current code and was not repeated here.
- *FY4B sets the grid.* Only hours with a strict-complete FY4B field are scored (13,471 of 26,304); winter is thinly represented.
- *Gauge resolution.* 0.5 mm tipping buckets cannot resolve 0.1–0.5 mm/h, which affects the 0.1 mm/h detection scores.

= Next steps

+ Re-run the benchmark without covariates, now that the configuration is on master, to settle whether they help.
+ Merge GPM and GSMaP over the full 2022–2024 record, without FY4B's gaps.
+ Produce gridded maps of the best fusion for the heavy-rain events already selected.

#pagebreak()

#heading(numbering: none)[Appendix — Reproduction]

#let prov(key) = { let hit = find(provenance, key: key); if hit == none { "?" } else { hit.value } }
Benchmark run, under `output/`:

#block(inset: (left: 1em))[
  #show raw: set text(size: 6.4pt)
  #raw(prov("run_directory"))
]

Commit #raw(prov("git_commit")) on branch #raw(prov("git_branch")). The run records uncommitted changes as #prov("git_dirty"): the only differences from the commit were untracked files that appeared in `Interim_results/` while the run was in progress (the heavy-rain spatial report and its figures); no tracked file differed, so the commit fully determines the code that ran. Every pooled metric is identical to the 2026-09-07 run of the same configuration and, for the twelve methods they share, to the 2026-09-11 run on master.

```sh
julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
    --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor
julia --project=. scripts/run_claim_reassessment.jl output/<run directory>
julia --project=. scripts/run_gauge_satellite_fusion.jl
py -3.13 scripts/plot_gauge_satellite_fusion.py
```

#table(
  columns: (46%, 54%),
  stroke: none,
  inset: (x: 5pt, y: 2.4pt),
  table.header(rule, [*File in* `output/gauge_satellite_fusion/`], [*Content*], thin),
  [`method_summary.csv`], [overall and detection skill per scheme, anchor and method],
  [`stratified_summary.csv`], [skill by rain intensity, distance to the nearest training gauge and year],
  [`fusion_vs_gauge_only.csv`], [every fusion configuration against ADW with bootstrap intervals],
  [`interpolation_ranking.csv`], [interpolators against each other],
  [`anchor_metrics.csv` \ `anchor_comparison.csv`], [anchors on identical cells, and against GPM],
  [`detection_skill_bootstrap.csv`], [POD, FAR and CSI against ADW with bootstrap intervals],
  [`merged_anchor_coefficients.csv` \ `blend_weights.csv`], [fitted merge weights and blending weights per fold],
  rule,
)
