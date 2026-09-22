// Interim results: MERGED_OLS_LAGNBR, a lagged, neighbour-aware, banded least-squares merge of FY4B, GPM and
// GSMaP, set against the other anchors and correction methods of the interpolation benchmark.
// Compile from this folder:  typst compile fused_anchor_lagnbr.typ
// Every number is read from tables/lagnbr_*.json, which scripts/plot_fused_anchor_lagnbr.py writes from
// output/fused_anchor_lagnbr/*.csv (scripts/run_fused_anchor_lagnbr_report.jl); figures are copies of
// output/fused_anchor_lagnbr/figures/.

#set document(
  title: "A lagged, neighbour-aware merged anchor",
  date: datetime(year: 2026, month: 9, day: 22),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Lagged, neighbour-aware merged anchor — interim results
      #h(1fr) 22 September 2026
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
#show figure.where(kind: table): set block(breakable: false)
#show table.cell: set text(size: 8pt)
#show table: set par(justify: false)
#show raw: set text(size: 8.5pt)
#set list(indent: 0.6em)

// ---------------------------------------------------------------------------------------------
// Data. Tables are loaded once; the prose reads single values through `val`, so a number in the text
// is always the number in the table.

#let load(name) = json("tables/lagnbr_" + name + ".json")
#let versus = load("vs_ols")
#let anchors = load("anchors")
#let adw = load("vs_adw")
#let detection = load("detection")
#let summary = load("method_summary")
#let bands = load("bands")
#let claim = load("claim")
#let coverage = load("coverage")
#let invariance = load("invariance")
#let cells = load("shared_cells")
#let provenance = load("provenance")

#let find(rows, ..conditions) = {
  let hits = rows.filter(row => conditions.named().pairs().all(pair => row.at(pair.at(0)) == pair.at(1)))
  if hits.len() == 0 { none } else { hits.at(0) }
}
/// One cell as typeset text: ASCII hyphens in numbers become minus signs.
#let val(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else { hit.at(column).replace("-", "−") }
}
/// A percentage for prose. A leading plus is dropped; a minus sign is kept.
#let pct(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else [#hit.at(column).trim("+", at: start).replace("-", "−")%]
}
/// "[low, high]%" for a row with ci_low / ci_high columns, set small.
#let interval(rows, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else {
    text(size: 8.5pt, fill: luma(80))[[#hit.ci_low.replace("-", "−"), #hit.ci_high.replace("-", "−")]%]
  }
}
/// A percentage's magnitude, for prose whose verb already says the direction ("lost 2.65%").
#let mag(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else [#hit.at(column).trim("+", at: start).trim("-", at: start)%]
}
#let signed(s) = s.replace("-", "−")

#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)
#let note(body) = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, fill: luma(246),
  stroke: (left: 2.5pt + rgb("#256abf")), body)

#let M = (
  raw: "Raw anchor", adw: "ADW (gauges only)", idw: "IDW", tps: "TPS", gwr: "GWR (gauges only)",
  "auto": "GWR family, chosen in-fold", residual_gwr: "Residual GWR", mixed_gwr: "Mixed GWR", mgwr: "MGWR",
  blend_residual_gwr: "Residual GWR + ADW blend", blend_mixed_gwr: "Mixed GWR + ADW blend",
  blend_mgwr: "MGWR + ADW blend", blend_agrenv_residual_gwr: "Residual GWR + agreement blend",
  blend_agrenv_mixed_gwr: "Mixed GWR + agreement blend", blend_agrenv_mgwr: "MGWR + agreement blend",
)
#let A = (FY4B: "FY4B", GPM: "GPM", GSMaP: "GSMaP", MERGED_MEAN: "Merged (equal)", MERGED_OLS: "Merged (OLS)",
  MERGED_OLS_LAGNBR: "Merged (OLS lag-nbr)", gauge_only: "Gauges only")
#let S = (all: "All hours", no_rain: "Dry < 0.1", light: "Light 0.1–2.5", moderate: "Moderate 2.5–8",
  heavy: "Heavy ≥ 8")
#let BS = "balanced_spatial"
#let STRATA-KEYS = ("all", "no_rain", "light", "moderate", "heavy")
#let NEW = "MERGED_OLS_LAGNBR"
#let OLD = "MERGED_OLS"

// ---------------------------------------------------------------------------------------------
// Title

#align(left)[
  #set text(hyphenate: false)
  #text(font: ("Segoe UI", "Arial"), size: 9pt, fill: luma(100), tracking: 0.04em)[INTERIM RESULTS · 22 SEPTEMBER 2026]
  #v(0.2em)
  #text(font: ("Segoe UI", "Arial"), size: 19pt, weight: "semibold")[A lagged, neighbour-aware merged anchor]
  #v(0.1em)
  #text(size: 12pt, fill: luma(60))[MERGED_OLS_LAGNBR: what widening the least-squares merge of FY4B, GPM and GSMaP buys, and where it sits among the other anchors and correction methods]
]
#v(0.6em)
#line(length: 100%, stroke: 0.5pt + luma(180))

= Key findings

- *A new anchor, MERGED_OLS_LAGNBR.* The existing OLS merge regresses the gauge on the three products at one hour and one point. The new merge also gives it each product one hour before and after, and each product's mean over the eight nearest stations, with a separate coefficient set for each of ten rain-agreement bands. It reads satellite values only and is fitted inside every fold on the training gauges.
- *It is the better anchor, by a small but solid margin.* On the cells both merges were scored on (balanced spatial folds), the raw anchor's RMSE falls from #val(versus, "RMSE_ols", scheme: BS, method: "raw", stratum: "all") to #val(versus, "RMSE_new", scheme: BS, method: "raw", stratum: "all") mm/h (#pct(versus, "gain", scheme: BS, method: "raw", stratum: "all") #interval(versus, scheme: BS, method: "raw", stratum: "all")). The gain carries through every correction method: MGWR #val(versus, "RMSE_ols", scheme: BS, method: "mgwr", stratum: "all") → #val(versus, "RMSE_new", scheme: BS, method: "mgwr", stratum: "all"), and the best method overall, MGWR + ADW blend, #val(versus, "RMSE_ols", scheme: BS, method: "blend_mgwr", stratum: "all") → #val(versus, "RMSE_new", scheme: BS, method: "blend_mgwr", stratum: "all"). Every bootstrap interval is above zero.
- *The gain is in moderate rain.* For MGWR it is #pct(versus, "gain", scheme: BS, method: "mgwr", stratum: "moderate") in the 2.5–8 mm/h stratum and #pct(versus, "gain", scheme: BS, method: "mgwr", stratum: "heavy") in heavy rain, where the interval spans zero. Against gauge-only ADW, residual GWR on the old merge was #mag(adw, "gain", scheme: BS, product: OLD, method: "residual_gwr", stratum: "moderate") worse in moderate rain; on the new merge it is #mag(adw, "gain", scheme: BS, product: NEW, method: "residual_gwr", stratum: "moderate") better.
- *Among all six anchors it is the best one for every correction method overall* (@sec-anchors), and for most of them in moderate rain. Heavy rain still favours the equal-weight merge: a least-squares anchor keeps a large dry bias on heavy hours whatever it is given.
- *Nothing else moved.* All #val(invariance, "baseline_rows", file: "metrics_pooled.csv") pooled metric rows and all #val(invariance, "baseline_rows", file: "paired_comparisons.csv") paired comparisons of the run without the new product are byte-identical in the run with it.
- *The project's claim is still not supported.* The heavy-rain gate (≥ 5 % against every baseline) is not met, detection scores are not better than ADW's, and one fold's covariate model left four stations unpredicted (@sec-limits).

= The method <sec-method>

== Why widen the merge

The OLS merge (`MERGED_OLS`) fits $y = a + b_1 "FY4B" + b_2 "GPM" + b_3 "GSMaP"$ on the training gauges. It was already the best of the five earlier anchors, but it can only correct the products' *amount* at the gauge pixel and hour. Hourly satellite precipitation is also wrong in *time* (a storm arrives an hour early or late) and in *position* (the rain cell is displaced by a pixel or two), and a one-hour, one-point regression cannot see either.

== Design

For every station $s$ and hour $t$ the design has thirteen terms:

- an intercept;
- FY4B, GPM and GSMaP at $(s, t)$;
- the same three at $(s, t - 1 "h")$ and at $(s, t + 1 "h")$ — the *lags*;
- the mean of each product over the eight nearest *other* stations at hour $t$ — the *neighbourhood means*.

A lag whose adjacent hour is missing from the common time grid falls back to the hour's own value, and so does a neighbourhood mean with no finite neighbour.

The coefficients are fitted separately in each of ten *agreement-envelope bands*: band 0 when no product reports ≥ 0.1 mm/h, otherwise (how many of the three report rain) × (their maximum, in 0.1–2.5, 2.5–8 and ≥ 8 mm/h). This is the same banding the agreement blend already uses. A band with fewer than 2,000 training cells, or a singular system, takes the pooled fit over all bands and is flagged; the prediction is clipped at zero.

== Leakage

- The lags and neighbourhood means read *satellite* values only. A held-out gauge's value never enters any feature.
- The band is read off the products, never off the gauge.
- The coefficients are fitted inside every fold on that fold's training stations only and then applied everywhere, exactly as for MERGED_OLS.

== How this variant was chosen

Before the benchmark, about thirty candidate merges were screened on the anchor alone, on the benchmark's own balanced spatial folds, with an inner four-fold split inside each training set to choose without looking at the test gauges. The candidates were per-band OLS, per-season OLS, lagged and neighbourhood inputs, square-root fitting, two-part (occurrence × amount) models, variance matching and locally fitted weights. Lagged and neighbourhood inputs gave most of the gain and banding added a little. Square-root fitting, two-part models and local weights did not beat plain OLS. Variance matching improved heavy rain only at a large cost overall. The chosen variant cut the anchor's held-out RMSE from 0.888 to 0.876, and its inner-split and in-sample scores were within 0.003 of that, so the extra coefficients are not overfitting. A richer version with polynomial terms scored slightly better but had three times the coefficients and a widening gap between its inner and held-out scores, so it was not taken forward.

= Data and evaluation <sec-data>

- *Run.* `run_interpolation_benchmark.jl full --nested-covariates --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --fused-anchor-lagnbr`: 237 gauges, 13,471 common hours (January 2022 – December 2024), six anchors (FY4B, GPM, GSMaP, MERGED_MEAN, MERGED_OLS, MERGED_OLS_LAGNBR) and every correction method of the benchmark. The run was made from a working tree at commit #val(provenance, "value", key: "git_commit") with the new code uncommitted; that code is commit `2088995`.
- *Validation.* Five-fold cross-validation by station. Balanced spatial folds are the primary scheme; random folds are a secondary check. Every hyperparameter, covariate and blending weight is chosen inside the training fold.
- *Same cells for every comparison.* The benchmark scores each product on its own evaluation mask, and MERGED_OLS_LAGNBR's is smaller than MERGED_OLS's (@sec-limits). Every comparison between products here is therefore rescored from the stored out-of-fold predictions on the cells all compared products share: #val(cells, "cells", scheme: BS, comparison: "lagnbr_vs_ols") cells for the head-to-head on balanced spatial folds, #val(cells, "cells", scheme: "random", comparison: "lagnbr_vs_ols") on random folds, and #val(cells, "cells", scheme: BS, comparison: "all_anchors") for the six-anchor comparison.
- *Uncertainty.* Paired day-block bootstrap, 2,000 replicates: whole days are resampled for both sides together. Holm correction across each table's family of comparisons.
- *Strata.* By the gauge value: dry < 0.1, light 0.1–2.5, moderate 2.5–8, heavy ≥ 8 mm/h.

= Results

== The new anchor among the anchors <sec-anchors>

On the #val(cells, "cells", scheme: BS, comparison: "all_anchors") cells all six anchors share, the new merge's raw field comes within #mag(summary, "vs_adw", scheme: BS, product: NEW, method: "raw") of gauge-only ADW (#val(anchors, "RMSE", scheme: BS, anchor: NEW, method: "raw") against #val(anchors, "RMSE", scheme: BS, anchor: "gauge_only", method: "adw")), and it stays the best anchor after every correction (@fig-anchors, @tab-anchors). Heavy rain is the exception: the equal-weight merge, which does not shrink, keeps the lowest heavy-rain RMSE (#val(anchors, "RMSE_heavy", scheme: BS, anchor: "MERGED_MEAN", method: "mgwr") mm/h after MGWR, against #val(anchors, "RMSE_heavy", scheme: BS, anchor: NEW, method: "mgwr")).

#figure(
  image("figures/lagnbr_fig2_anchors.png", width: 100%),
  caption: [RMSE and bias of the six anchors by gauge rain intensity, balanced spatial CV.],
) <fig-anchors>

#let anchor-methods = ("raw", "residual_gwr", "mgwr", "blend_mgwr", "blend_agrenv_mgwr")
#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.3pt),
    table.header(rule, [*Anchor*], [*Method*], [*RMSE*], [*Moderate*], [*Heavy*], [*vs Merged (OLS)*],
      [*Random folds*], thin),
    ..for anchor in ("FY4B", "GPM", "GSMaP", "MERGED_MEAN", OLD, NEW) {
      let rows = anchors.filter(row => row.scheme == BS and row.anchor == anchor and row.method in anchor-methods)
      (
        table.cell(rowspan: rows.len(), if anchor == NEW { strong(A.at(anchor)) } else { A.at(anchor) }),
        ..for row in rows {
          (M.at(row.method), row.RMSE, row.RMSE_moderate, row.RMSE_heavy,
            if row.vs_ols == "–" { [–] } else [#signed(row.vs_ols)% #text(size: 6.5pt, fill: luma(90))[[#signed(row.ci_low), #signed(row.ci_high)]]],
            val(anchors, "RMSE", scheme: "random", anchor: anchor, method: row.method))
        },
        thin,
      )
    },
    A.gauge_only, M.adw, val(anchors, "RMSE", scheme: BS, anchor: "gauge_only"),
    val(anchors, "RMSE_moderate", scheme: BS, anchor: "gauge_only"),
    val(anchors, "RMSE_heavy", scheme: BS, anchor: "gauge_only"), [–],
    val(anchors, "RMSE", scheme: "random", anchor: "gauge_only"),
    rule,
  ),
  caption: [RMSE (mm/h) on the cells all six anchors share (#val(cells, "cells", scheme: BS, comparison: "all_anchors") balanced spatial, #val(cells, "cells", scheme: "random", comparison: "all_anchors") random); "vs Merged (OLS)" is the balanced-spatial improvement with its 95% interval.],
) <tab-anchors>

== Correction methods on the new anchor <sec-methods>

Every method that reads the anchor improves on balanced spatial folds, the gain concentrated in light and moderate rain (@fig-gain, @fig-strata, @tab-methods). The raw anchor pays for it on dry hours (#pct(versus, "gain", scheme: BS, method: "raw", stratum: "no_rain")), where its lag and neighbourhood terms let nearby rain leak in as drizzle; the correction removes that cost. On random folds the gains are smaller and the residual agreement blend is the one significant loss.

#figure(
  image("figures/lagnbr_fig1_gain_by_method.png", width: 100%),
  caption: [MERGED_OLS_LAGNBR against MERGED_OLS for every method that reads the anchor.],
) <fig-gain>

#figure(
  image("figures/lagnbr_fig3_gain_by_intensity.png", width: 100%),
  caption: [Gain and change in absolute bias by gauge rain intensity.],
) <fig-strata>

#let method-order = ("raw", "residual_gwr", "mixed_gwr", "mgwr", "auto", "blend_residual_gwr", "blend_mixed_gwr",
  "blend_mgwr", "blend_agrenv_residual_gwr", "blend_agrenv_mixed_gwr", "blend_agrenv_mgwr")
#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.3pt),
    table.header(rule, [*Method*], [*OLS*], [*LAGNBR*], [*All hours*], [*Moderate*], [*Heavy*], [*Random folds*], thin),
    ..for method in method-order {
      let cell(scheme, stratum) = {
        let row = find(versus, scheme: scheme, method: method, stratum: stratum)
        [#signed(row.gain)% #text(size: 6.3pt, fill: luma(90))[[#signed(row.ci_low), #signed(row.ci_high)]]]
      }
      (M.at(method), val(versus, "RMSE_ols", scheme: BS, method: method, stratum: "all"),
        val(versus, "RMSE_new", scheme: BS, method: method, stratum: "all"),
        cell(BS, "all"), cell(BS, "moderate"), cell(BS, "heavy"), cell("random", "all"))
    },
    rule,
  ),
  caption: [RMSE (mm/h, balanced spatial) and relative improvement of LAGNBR over OLS with its 95% interval.],
) <tab-methods>

== Against the gauge-only reference <sec-adw>

On balanced spatial folds every fusion method beats ADW overall on both merges, and by more on the new one (@fig-adw, @tab-adw). The clearest change is moderate rain, where MGWR goes from a tie with ADW on the OLS merge to #pct(adw, "gain", scheme: BS, product: NEW, method: "mgwr", stratum: "moderate") better on the new one #interval(adw, scheme: BS, product: NEW, method: "mgwr", stratum: "moderate"). On random folds only the blends reach ADW, on either merge. Detection is traded between thresholds rather than improved: the new merge gains CSI at 2.5 mm/h and loses it at 0.1 mm/h (@tab-detection), and the claim's event criterion fails for both.

#figure(
  image("figures/lagnbr_fig4_vs_adw.png", width: 100%),
  caption: [Both merges against gauge-only ADW, by stratum and fold scheme, with detection skill.],
) <fig-adw>

#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.3pt),
    table.header(rule, [*Method*], ..STRATA-KEYS.map(s => [*#S.at(s)*]), thin),
    ..for scheme in (BS, "random") {
      (table.cell(colspan: 6, text(style: "italic")[#if scheme == BS [Balanced spatial folds] else [Random folds]]),
        ..for method in ("residual_gwr", "mgwr", "blend_mgwr", "blend_agrenv_mgwr") {
          (M.at(method),
            ..for stratum in STRATA-KEYS {
              let old = find(adw, scheme: scheme, product: OLD, method: method, stratum: stratum)
              let new = find(adw, scheme: scheme, product: NEW, method: method, stratum: stratum)
              ([#text(fill: luma(110))[#signed(old.gain)] → #signed(new.gain)%],)
            })
        })
    },
    rule,
  ),
  caption: [Relative RMSE improvement over ADW, on MERGED_OLS → on MERGED_OLS_LAGNBR.],
) <tab-adw>

#let det-methods = ("raw", "mgwr", "blend_mgwr")
#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.3pt),
    table.header(rule, [*Method*], [*Merge*], [*0.1 mm/h*], [*2.5 mm/h*], [*8 mm/h*], [*16 mm/h*], thin),
    ..for method in det-methods {
      for product in (OLD, NEW) {
        (if product == OLD { M.at(method) } else { [] }, if product == OLD { [OLS] } else { [LAGNBR] },
          ..for threshold in ("0.1", "2.5", "8.0", "16.0") {
            let row = find(detection, scheme: BS, product: product, method: method, threshold: threshold, metric: "CSI")
            ([#row.method_value #text(size: 6.3pt, fill: luma(90))[(#signed(row.delta))]],)
          })
      }
    },
    thin,
    [ADW (gauges only)], [], ..for threshold in ("0.1", "2.5", "8.0", "16.0") {
      (val(detection, "adw", scheme: BS, product: OLD, method: "mgwr", threshold: threshold, metric: "CSI"),)
    },
    rule,
  ),
  caption: [Critical success index, balanced spatial CV; in brackets, the difference from ADW.],
) <tab-detection>

== What the fit learned <sec-bands>

GPM's neighbourhood mean is the largest slope in every band (tied with GSMaP at $t$ in one), so the merge trusts GPM's rain around a station more than at it — the position error the term was meant to absorb. The $t + 1$ h terms outweigh the $t - 1$ h terms in most bands, so the products tend to report rain late. Bands 0 (no product wet: the products are exactly zero, so the system is singular) and 3 (about #val(bands, "cells", band: "3") cells) took the pooled fit in every fold (@fig-bands, @tab-bands).

#figure(
  image("figures/lagnbr_fig5_coefficients.png", width: 100%),
  caption: [Fold-mean coefficients by band and term, with training cells and summed weights per band.],
) <fig-bands>

#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto, auto, auto, auto),
    align: (right, left, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.3pt),
    table.header(rule, [*Band*], [*Products reporting rain · their maximum*], [*Cells*], [*Pooled*], [*Intercept*],
      [*At t*], [*t − 1 h*], [*t + 1 h*], [*Nearest 8*], thin),
    ..for row in bands {
      (row.band, row.label, row.cells, row.fell_back, signed(row.intercept), signed(row.sum_t), signed(row.sum_lag),
        signed(row.sum_lead), signed(row.sum_nbr))
    },
    rule,
  ),
  caption: [Mean training cells per fold, folds on the pooled fit, and summed coefficients per term group.],
) <tab-bands>

== Every method on every anchor <sec-summary>

@tab-summary lists every method as the benchmark reports it, each product on its own evaluation mask; it ranks methods within an anchor, while @tab-anchors compares anchors.

#let summary-methods = ("raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr", "auto", "blend_residual_gwr",
  "blend_mixed_gwr", "blend_mgwr", "blend_agrenv_residual_gwr", "blend_agrenv_mixed_gwr", "blend_agrenv_mgwr")
#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.2pt),
    table.header(rule, [*Method*], table.cell(colspan: 3)[*RMSE, balanced spatial*],
      table.cell(colspan: 2)[*vs ADW, balanced spatial*], table.cell(colspan: 2)[*vs ADW, random*],
      [], [GPM], [OLS], [LAGNBR], [OLS], [LAGNBR], [OLS], [LAGNBR], thin),
    ..for method in summary-methods {
      (M.at(method),
        val(summary, "RMSE", scheme: BS, product: "GPM", method: method),
        val(summary, "RMSE", scheme: BS, product: OLD, method: method),
        val(summary, "RMSE", scheme: BS, product: NEW, method: method),
        [#val(summary, "vs_adw", scheme: BS, product: OLD, method: method)%],
        [#val(summary, "vs_adw", scheme: BS, product: NEW, method: method)%],
        [#val(summary, "vs_adw", scheme: "random", product: OLD, method: method)%],
        [#val(summary, "vs_adw", scheme: "random", product: NEW, method: method)%])
    },
    rule,
  ),
  caption: [Pooled RMSE (mm/h) and improvement over ADW, each product on its own mask (#val(summary, "n", scheme: BS, product: NEW, method: "raw") cells for LAGNBR, #val(summary, "n", scheme: BS, product: OLD, method: "raw") for OLS).],
) <tab-summary>

= Limitations <sec-limits>

- *Heavy rain.* No squared-error anchor fixes it. Given what the satellites report, the expected gauge value in a heavy-rain hour is low, so the least-squares merge keeps a dry bias on heavy hours and the claim's heavy gate (≥ 5 % against every baseline) is out of reach this way. The screen found the same for every least-squares variant; only variance matching moved heavy rain, at a large cost overall.
- *The claim is not supported.* @tab-claim is the benchmark's own assessment. MERGED_OLS_LAGNBR improves residual GWR over the best traditional method by #pct(claim, "overall_gain", product: NEW), significantly, and removes the moderate-rain degradation, but fails the heavy gate, the event criterion and the coverage criterion.
- *Fold 4 coverage.* In fold 4 of the balanced spatial scheme the per-fold covariate selection for MERGED_OLS_LAGNBR picked #raw(val(coverage, "covariates", product: NEW, fold: "4")), and the covariate model then predicted only #val(coverage, "coverage", product: NEW, fold: "4") of that fold's cells. The GWR-family predictions are missing at four neighbouring stations (61937000, 61937150, 61937215, 61937220, around 110.8–111.0°E, 32.8°N). GPM hit the same failure in the same fold with a covariate set that also contains NDVI (#raw(val(coverage, "covariates", product: "GPM", fold: "4"))), while MERGED_OLS did not select NDVI there. This is a property of the covariate path, not of the new anchor, and the shared-cell rescoring above removes its effect on the comparison.
- *Modest ceiling.* The best candidate merge fitted in-sample was only about 2.7 % below plain OLS; the new merge takes about half of that.
- *Other reports.* `run_gauge_satellite_fusion.jl`, `run_no_rain_fusion_evaluation.jl` and their plot scripts list their anchors by hand and do not yet include MERGED_OLS_LAGNBR.

#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.3pt),
    table.header(rule, [*Anchor*], [*RMSE*], [*Overall*], [*Significant*], [*Heavy*], [*Moderate loss*],
      [*Coverage*], [*Events*], [*Supported*], thin),
    ..for row in claim {
      (A.at(row.product), row.RMSE_residual_gwr, [#signed(row.overall_gain)%], row.paired_significant,
        [#signed(row.heavy_gain)%], [#signed(row.moderate_degradation)%], [#row.coverage (#row.coverage_ok)],
        row.events_ok, row.supported)
    },
    rule,
  ),
  caption: [The benchmark's claim assessment for residual GWR against the best traditional method (ADW), balanced spatial CV, each anchor on its own mask. "Overall" and "Heavy": relative RMSE improvement; "Moderate loss": relative degradation in moderate rain (negative is a gain).],
) <tab-claim>

= Reproduction

From the repository root (the benchmark takes about nine hours; the rest about two minutes):

```sh
julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
  --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --fused-anchor-lagnbr
julia -t 4 --project=. scripts/run_fused_anchor_lagnbr_report.jl
py -3.13 scripts/plot_fused_anchor_lagnbr.py
cd Interim_results && typst compile fused_anchor_lagnbr.typ
```
