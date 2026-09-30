// Paper plan with draft figures: gauge, GPM and GSMaP precipitation in a sparse-gauge mountain area.
// Compile from the repository root:  typst compile paper/paper_plan.typ
// Figures are copies of the `--gpm-gsmap` outputs of scripts/plot_heavy_rain_events.py,
// scripts/plot_satellite_temporal_evaluation.py and scripts/plot_gauge_satellite_fusion.py.
// Numbers are read from the CSVs beside those figures in output/.

#set document(
  title: "Paper plan: satellite precipitation and gauge fusion (GPM, GSMaP)",
  date: datetime(year: 2026, month: 9, day: 27),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Paper plan with draft figures — GPM and GSMaP
      #h(1fr) 27 September 2026
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
#show table.cell: set text(size: 8.5pt)
#show table: set par(justify: false)
#set list(indent: 0.6em)

#let note(body) = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, fill: luma(246),
  stroke: (left: 2.5pt + rgb("#2a78d6")), body)
#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)

#align(center)[
  #text(size: 17pt, weight: "semibold", font: ("Segoe UI", "Arial"))[Paper plan with draft figures]
  #v(0.2em)
  #text(size: 11.5pt, fill: luma(60))[Can GPM and GSMaP stand in for gauges where the network is sparse, and what does fusing the two gain?]
]

#note[
  *Scope.* This plan covers two satellite products, GPM IMERG and GSMaP. FY4B is left out of the paper.
  It answers the three questions set for the study and the two methodological questions, in eight
  figures. Five of them are shown below as working drafts built from finished analyses; the other
  three are summarised in the text. Final styling, panel letters and captions come later.
]

= Aim, data and storyline

In a sparse gauge network, satellites give full spatial coverage but with errors. The paper asks
three questions in order. (Q1) Does the satellite field reproduce *where* heavy rain falls? (Q2)
Does it reproduce *how rain evolves in time* at each intensity? (Q3) What does fusing gauges with
the satellite gain, and does it fix the moderate and heavy rain that flood forecasting depends on?
Two methodological answers follow: which interpolation and fusion methods are best, and which
product (or a merge of both) should be the anchor.

*Data.* 237 hourly gauges in the Shiyan study area (data quantised to 0.5 mm), and GPM IMERG and
GSMaP hourly at 0.1°, sampled at the gauge pixels, for 2022–2024. Days run 08–08 Beijing time. Local
relief from the Copernicus GLO-30 DEM defines the landform regions.

*Figures.* Eight are planned. Five are shown below as drafts (Figs. 2, 3, 4, 6, 7); the other three
are summarised in the text.
+ Study area, gauges and landform regions.
+ Six heavy-rain events: gauge, GPM and GSMaP maps (Q1).
+ Event-scale spatial scores (Q1).
+ Hourly skill by intensity class, by season and region (Q2).
+ Rain events: detection, timing and volume by class (Q2).
+ Gauge-only interpolators; fusion method × anchor (Q3, best method).
+ Where fusion helps: intensity, distance to gauges, detection (Q3).
+ Which product, and does merging help (best product).

*Study area (Fig. 1).* The 237 gauges sit in hills and mountains; none is on a plain. By local relief,
129 gauges are below 500 m, 97 are at 500–1000 m and 11 are at 1000 m or more (too few to interpret
alone). Fig. 1 will be a map with relief, gauges, a location inset and river basins.

= Q1 — Spatial heterogeneity in heavy-rain events

*Design.* Heavy-rain days are 08–08 days on which more than three gauges exceed 50 mm. Of 56 such
days, six are selected, one per rain system: four widespread and two localized. At each gauge, the
satellite 24 h total is compared with the gauge total, using spatial correlation, spatial CV ratio,
rain-centre shift, mean bias, RMSE and CSI at 50 mm.

#figure(
  image("figures/fig02_event_maps.png", width: 100%),
  caption: [Draft Fig. 2. The six heavy-rain events (columns). Top row: gauge 24 h totals (08–08 BJT).
  Lower rows: satellite minus gauge at each gauge, from the pixel sampled there; red means the satellite is
  drier, blue wetter, grey within ±10 mm. MB is the mean bias (mean of satellite − gauge over the
  event's gauges).],
)

#figure(
  image("figures/fig03_event_metrics.png", width: 100%),
  caption: [Draft Fig. 3. Spatial scores per event and product.],
)

*Key findings.*
- Both products place heavy rain about equally well. Mean spatial r over the six events is 0.55
  for both GPM and GSMaP, and 0.6–0.8 on four events. Both fail on 2022-06-26 (r ≤ 0.12).
- Amounts differ sharply. GPM is nearly unbiased (MB −2.0 mm, RMSE 20.5 mm). GSMaP
  overestimates (MB +18.9 mm, RMSE 39.9 mm), up to +73 mm on 2023-08-26.
- Localized storms are underestimated by both products: CSI at 50 mm is 0 except for GPM on
  2024-07-30. The 0.1° pixel smooths out small intense cores.

= Q2 — Temporal evolution of rainfall

*Design.* Only gauge-wet hours are scored; no-rain hours are excluded. The hourly intensity classes
are 0.1–2, 2–4, 4–8, 8–20 and ≥ 20 mm h⁻¹. Results are stratified by season (MAM, JJA, SON, DJF)
and by landform region, with day-block bootstrap intervals.

#figure(
  image("figures/fig04a_intensity_season.png", width: 100%),
  caption: [Draft Fig. 4. Hourly skill by gauge intensity class and season. The paper version adds the
  landform-region split.],
)

*Key findings.*
- Underestimation grows with intensity for both products. GPM relative bias runs from −11% (light)
  to −38% (4–8), −68% (8–20) and −85% (≥ 20 mm h⁻¹). GSMaP runs from +16% to +8%, −17%, −59% and −80%.
  At ≥ 8 mm h⁻¹, 83–99% of hours fall into a lower class.
- Detection is not the problem (Fig. 5, rain events separated by at least 3 dry hours). Rain events peaking at 2–4 mm h⁻¹ or more are detected 88–96% of the time, and
  rain-centre timing errors stay within ±0.5 h. The error is in peak *magnitude*.
- Winter (DJF) is where both fail: bias about −80% and only 14–18% of wet hours detected. Summer
  light rain is overestimated (GPM +39%, GSMaP +52%).
- Landform regions differ little. Between relief below 500 m and 500–1000 m, the class biases differ by
  at most 8 percentage points, far less than between classes or seasons. The ≥ 1000 m region has 11 gauges
  and is not interpreted.
- Areal-mean series: hourly r is 0.83 for GPM and 0.66 for GSMaP (daily: 0.91 and 0.80). The overall
  bias is +6% for GPM and +33% for GSMaP.

= Q3 — Benefits of gauge–satellite fusion

*Design.* Every method is scored at gauges held out of the fit, using balanced spatial 5-fold
cross-validation over 13,471 hours. Gauge-only interpolators are IDW, ADW, TPS and GWR. Fusion methods use the
satellite as an anchor (residual GWR, mixed GWR, MGWR), optionally blended with ADW where the
satellite reports rain.

#figure(
  image("figures/fig06b_fusion_methods.png", width: 100%),
  caption: [Draft Fig. 6. RMSE improvement of each fusion method over gauge-only ADW, per anchor, for
  all hours and heavy hours. The paper version adds a panel comparing the gauge-only interpolators.],
)

#figure(
  image("figures/fig07_where_fusion_helps.png", width: 100%),
  caption: [Draft Fig. 7. Where the satellite adds to the gauges: by rain intensity, by distance to
  the nearest training gauge, and in detection skill (CSI).],
)

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    rule,
    table.header([*Anchor*], [*Method*], [*RMSE*], [*r*], [*CSI 2.5*], [*vs raw*], [*vs ADW*]),
    thin,
    [—], [ADW (gauges only)], [0.871], [0.486], [0.280], [—], [0.0%],
    [—], [IDW (gauges only)], [0.874], [0.483], [0.279], [—], [−0.3%],
    [—], [TPS (gauges only)], [0.931], [0.437], [0.281], [—], [−6.9%],
    [—], [GWR (gauges only)], [0.922], [0.422], [0.250], [—], [−5.9%],
    thin,
    [GPM], [raw], [0.962], [0.422], [0.253], [0.0%], [−10.4%],
    [GPM], [MGWR + agreement blend], [0.844], [0.528], [0.317], [+12.3%], [+3.1%],
    [GSMaP], [raw], [1.219], [0.369], [0.233], [0.0%], [−40.0%],
    [GSMaP], [MGWR + agreement blend], [0.853], [0.520], [0.313], [+30.1%], [+2.1%],
    [Merged (OLS)], [raw], [0.888], [0.453], [0.220], [0.0%], [−1.4%],
    [Merged (OLS)], [MGWR + agreement blend], [0.844], [0.532], [0.308], [+5.0%], [+3.6%],
    rule,
  ),
  caption: [Draft headline table. Held-out gauges, balanced spatial CV; RMSE in
  mm h⁻¹; improvements are RMSE reductions. Gauge-only rows are on GPM's evaluation mask.],
)

*Key findings.*
- Fusion beats both inputs. With MGWR plus the agreement blend, RMSE falls by 12% from raw GPM and by
  30% from raw GSMaP, and it beats gauge-only ADW by 2–4%. Every blended variant on GPM and on the
  merged anchor beats ADW significantly, and so do the MGWR-based blends on GSMaP.
- The gain is concentrated where it matters for floods. On heavy hours (≥ 8 mm h⁻¹), MGWR alone
  improves on ADW by 6.0% (GPM) and 7.4% (GSMaP). Unblended regression, however, *hurts* light
  and dry hours, which is why the blend is needed.
- The gain grows with distance from the network. On GPM, the improvement over ADW is 1% within 20 km
  of a training gauge, 3% at 20–50 km and 7% at 50–100 km (GSMaP: 1%, 2%, 4%). This is the sparse-network case the study motivates.

= Methodological answers

Both answers rest on the headline table above. Fig. 8 in the paper will show raw vs fused skill per
anchor and the merge weights.

- *M1, best interpolator.* ADW, with IDW indistinguishable from it (0.871 vs 0.874 mm h⁻¹). TPS and
  gauge-only GWR are 6–7% worse. *Best fusion:* MGWR followed by a blend with ADW where the satellite
  reports rain. MGWR alone is best for heavy hours only.
- *M2, best product.* GPM. It has the lower raw error (RMSE 0.96 vs 1.22 mm h⁻¹), near-zero event bias
  and a better areal-mean timeline. GSMaP's large wet bias is largely removed by fusion, and once fused
  the two anchors differ by only about 1%.
- *Merging.* The OLS merge is a better raw field than either product (RMSE 0.888). After fusion it adds
  almost nothing over GPM (RMSE 0.844 for both). Merging is worth reporting, but it is not the main source of gain.

= Flood-forecasting extension (outlook)

Q2 shows that raw satellite rain loses about 60–85% of the volume at 8 mm h⁻¹ and above. That is the input
flood models are most sensitive to, so raw products should not drive them. Fusion recovers part of
this (Figs. 6 and 7), but underestimation at the highest intensities remains. The planned extension:

+ Aggregate raw, gauge-only and fused hourly fields to basin-mean rainfall for the gauged
  catchments during the six Q1 events and the wettest seasons.
+ Correlate basin rainfall (and antecedent totals) with observed hydrological-station discharge, and
  compare raw vs fused.
+ Only if the correlation improves clearly, drive a simple rainfall–runoff model as a demonstration.

*Missing input:* discharge data are not yet available in the project and must be obtained first.

= Open items before drafting

- Fig. 1: draw the study-area map (relief, gauges, location inset, basins).
- Fig. 7 and Table 3 use fusion intensity strata (0.1–2.5, 2.5–8, ≥ 8 mm h⁻¹) that differ from the Q2
  classes. Re-stratify the fusion results on the Q2 classes for consistency.
- The merged anchor still carries a small FY4B term (OLS weight 0.03, against 0.40 for GPM and 0.16
  for GSMaP). It stays for now; refit as a GPM + GSMaP merge before submission.
- Decide whether the newer stacked combination of MGWR and ADW (`stack_mgwr`) replaces the agreement
  blend as the headline fusion method.
- Diurnal-cycle results exist and can go to supplementary material.
