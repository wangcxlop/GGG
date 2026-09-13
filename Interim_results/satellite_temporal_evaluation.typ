// Interim results: temporal evaluation of FY4B, GPM and GSMaP against gauges, 2022-2024.
// Compile from this folder:  typst compile satellite_temporal_evaluation.typ
// Every number in the tables is generated from output/satellite_temporal_evaluation/*.csv
// (commit 60012c4, scripts/run_satellite_temporal_evaluation.jl); figures are copies of
// output/satellite_temporal_evaluation/figures/.

#set document(
  title: "Temporal evolution of rainfall in FY4B, GPM and GSMaP: interim results",
  date: datetime(year: 2026, month: 9, day: 13),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Temporal evaluation of satellite precipitation — interim results
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
// Long tables continue on the next page; the short ones are wrapped in `block(breakable: false)`.
#show figure.where(kind: table): set block(breakable: true)
#show table.cell: set text(size: 8pt)
#show table: set par(justify: false)
#show raw: set text(size: 8.5pt)
#set list(indent: 0.6em)

// ---------------------------------------------------------------------------------------------
// Tables. The row data below was generated from the CSVs; do not edit numbers by hand.

#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)

/// Two sample groups side by side: the hours FY4B covers (all three products on identical cells) and
/// the full record for GPM and GSMaP.
#let split-table(first, rows, first-width: 30%) = table(
  columns: (first-width, 1fr, 1fr, 1fr, 1fr, 1fr),
  align: (left + horizon, right + horizon, right + horizon, right + horizon, right + horizon, right + horizon),
  stroke: none,
  inset: (x: 4pt, y: 2.6pt),
  table.vline(x: 4, stroke: 0.4pt + luma(170)),
  table.header(
    rule,
    table.cell(rowspan: 2, align: left + bottom)[*#first*],
    table.cell(colspan: 3, align: center)[*Hours FY4B covers* (three products)],
    table.cell(colspan: 2, align: center)[*Full 2022–2024 record*],
    [*FY4B*], [*GPM*], [*GSMaP*], [*GPM*], [*GSMaP*],
    thin,
  ),
  ..rows,
  rule,
)

#let coverage_rows = (
  [All seasons], [13,243], [614], [214,803], [7.0], [26,304], [1,096], [381,393], [6.3],
  [MAM], [2,589], [127], [48,412], [8.0], [6,624], [276], [106,773], [7.0],
  [JJA], [6,091], [276], [95,937], [6.7], [6,624], [276], [104,491], [6.7],
  [SON], [3,160], [149], [62,770], [8.6], [6,552], [273], [126,988], [8.4],
  [DJF], [1,403], [62], [7,684], [2.4], [6,504], [271], [43,141], [3.0],
)

#let landform_rows = (
  table.cell(colspan: 5, fill: luma(235))[*Basic landform class*],
  [Plain (< 30 m)], [0], [–], [–], [–],
  [Hills (30–200 m)], [13], [142], [9.5], [274],
  [Small-relief mountain (200–500 m)], [116], [363], [20.4], [459],
  [Medium-relief mountain (500–1000 m)], [97], [674], [27.6], [813],
  [Large-relief mountain (≥ 1000 m)], [11], [1079], [33.5], [1081],
  table.cell(colspan: 5, fill: luma(235))[*Evaluation region*],
  [Relief < 500 m (hills + small-relief mountains)], [129], [341], [19.3], [440],
  [Relief 500–1000 m], [97], [674], [27.6], [813],
  [#text(fill: luma(110))[Relief ≥ 1000 m#super[a]]], [#text(fill: luma(110))[11]], [#text(fill: luma(110))[1079]], [#text(fill: luma(110))[33.5]], [#text(fill: luma(110))[1081]],
)

#let intensity_rows = (
  table.cell(colspan: 6, fill: luma(235))[*Gauge-wet station-hours (n)*],
  [Light (0.1–2)], [152,542], [152,542], [152,542], [286,433], [286,433],
  [Moderate (2–4)], [37,396], [37,396], [37,396], [60,202], [60,202],
  [Heavy (4–8)], [17,138], [17,138], [17,138], [25,266], [25,266],
  [Rainstorm (8–20)], [6,302], [6,302], [6,302], [7,866], [7,866],
  [Severe (≥ 20)], [1,425], [1,425], [1,425], [1,626], [1,626],
  [All wet hours], [214,803], [214,803], [214,803], [381,393], [381,393],
  table.cell(colspan: 6, fill: luma(235))[*Relative bias of the hourly amount (%) with 95% CI*],
  [Light (0.1–2)], [−17 #text(size: 6.5pt, fill: luma(90))[(−40, 12)]], [9 #text(size: 6.5pt, fill: luma(90))[(−2, 21)]], [33 #text(size: 6.5pt, fill: luma(90))[(14, 54)]], [−11 #text(size: 6.5pt, fill: luma(90))[(−19, −3)]], [16 #text(size: 6.5pt, fill: luma(90))[(3, 32)]],
  [Moderate (2–4)], [−44 #text(size: 6.5pt, fill: luma(90))[(−60, −26)]], [−15 #text(size: 6.5pt, fill: luma(90))[(−25, −6)]], [9 #text(size: 6.5pt, fill: luma(90))[(−15, 34)]], [−20 #text(size: 6.5pt, fill: luma(90))[(−27, −12)]], [8 #text(size: 6.5pt, fill: luma(90))[(−10, 26)]],
  [Heavy (4–8)], [−56 #text(size: 6.5pt, fill: luma(90))[(−68, −44)]], [−41 #text(size: 6.5pt, fill: luma(90))[(−49, −35)]], [−18 #text(size: 6.5pt, fill: luma(90))[(−41, 3)]], [−38 #text(size: 6.5pt, fill: luma(90))[(−44, −32)]], [−17 #text(size: 6.5pt, fill: luma(90))[(−34, −1)]],
  [Rainstorm (8–20)], [−74 #text(size: 6.5pt, fill: luma(90))[(−83, −64)]], [−69 #text(size: 6.5pt, fill: luma(90))[(−74, −64)]], [−59 #text(size: 6.5pt, fill: luma(90))[(−71, −47)]], [−68 #text(size: 6.5pt, fill: luma(90))[(−72, −64)]], [−59 #text(size: 6.5pt, fill: luma(90))[(−69, −49)]],
  [Severe (≥ 20)], [−84 #text(size: 6.5pt, fill: luma(90))[(−90, −79)]], [−86 #text(size: 6.5pt, fill: luma(90))[(−88, −84)]], [−81 #text(size: 6.5pt, fill: luma(90))[(−87, −74)]], [−85 #text(size: 6.5pt, fill: luma(90))[(−87, −83)]], [−80 #text(size: 6.5pt, fill: luma(90))[(−87, −73)]],
  [All wet hours], [−48 #text(size: 6.5pt, fill: luma(90))[(−61, −34)]], [−30 #text(size: 6.5pt, fill: luma(90))[(−37, −25)]], [−11 #text(size: 6.5pt, fill: luma(90))[(−27, 5)]], [−32 #text(size: 6.5pt, fill: luma(90))[(−37, −28)]], [−11 #text(size: 6.5pt, fill: luma(90))[(−22, 1)]],
  table.cell(colspan: 6, fill: luma(235))[*Detected as rain, satellite ≥ 0.1 mm (%)*],
  [Light (0.1–2)], [27], [58], [58], [51], [52],
  [Moderate (2–4)], [49], [82], [82], [81], [81],
  [Heavy (4–8)], [59], [85], [86], [86], [87],
  [Rainstorm (8–20)], [62], [83], [85], [83], [85],
  [Severe (≥ 20)], [70], [85], [89], [85], [89],
  [All wet hours], [34], [65], [66], [59], [60],
  table.cell(colspan: 6, fill: luma(235))[*Satellite in the same intensity class (%)*],
  [Light (0.1–2)], [18], [45], [44], [41], [40],
  [Moderate (2–4)], [10], [21], [17], [20], [16],
  [Heavy (4–8)], [9], [21], [18], [21], [20],
  [Rainstorm (8–20)], [10], [12], [12], [13], [12],
  [Severe (≥ 20)], [4], [1], [6], [1], [6],
  [All wet hours], [16], [38], [36], [36], [34],
  table.cell(colspan: 6, fill: luma(235))[*Satellite in a lower class, including no rain (%)*],
  [Light (0.1–2)], [73], [42], [42], [49], [48],
  [Moderate (2–4)], [81], [62], [62], [64], [62],
  [Heavy (4–8)], [82], [71], [67], [69], [65],
  [Rainstorm (8–20)], [88], [87], [83], [87], [83],
  [Severe (≥ 20)], [96], [99], [94], [99], [94],
  [All wet hours], [76], [50], [49], [54], [52],
)

#let false_alarm_rows = (
  table.cell(colspan: 6, fill: luma(235))[*Satellite ≥ 0.1 mm on gauge-dry hours (%)*],
  [All seasons], [6.1], [7.0], [7.2], [5.6], [6.3],
  [MAM], [6.3], [6.8], [6.8], [5.7], [6.7],
  [JJA], [6.5], [8.4], [7.6], [8.3], [7.5],
  [SON], [7.8], [6.9], [7.5], [6.3], [7.3],
  [DJF], [0.3], [1.5], [5.3], [1.9], [3.6],
  table.cell(colspan: 6, fill: luma(235))[*Satellite ≥ 0.5 mm (the gauge resolution) on gauge-dry hours (%)*],
  [All seasons], [3.9], [3.1], [3.1], [2.4], [2.6],
  [MAM], [3.8], [2.8], [3.0], [2.4], [3.1],
  [JJA], [4.1], [4.0], [3.3], [3.9], [3.3],
  [SON], [5.2], [2.6], [3.3], [2.4], [3.1],
  [DJF], [0.1], [0.6], [1.4], [0.8], [1.0],
  table.cell(colspan: 6, fill: luma(235))[*False-alarm ratio: gauge-dry share of satellite-wet hours (%)*],
  [All seasons], [70], [59], [59], [58], [61],
  [MAM], [76], [55], [55], [56], [59],
  [JJA], [70], [62], [61], [62], [61],
  [SON], [67], [53], [54], [52], [55],
  [DJF], [73], [71], [85], [82], [86],
  table.cell(colspan: 6, fill: luma(235))[*Share of satellite rain volume on gauge-dry hours (%)*],
  [All seasons], [60], [36], [31], [36], [33],
  [MAM], [74], [34], [36], [34], [35],
  [JJA], [62], [38], [30], [38], [30],
  [SON], [53], [30], [28], [30], [29],
  [DJF], [58], [61], [76], [76], [80],
)

#let event_rows = (
  table.cell(colspan: 6, fill: luma(235))[*Scored events (n)*],
  [All events], [23,955], [23,955], [23,955], [86,508], [86,508],
  [Light], [16,883], [16,883], [16,883], [57,195], [57,195],
  [Moderate], [3,376], [3,376], [3,376], [13,573], [13,573],
  [Heavy], [1,989], [1,989], [1,989], [9,021], [9,021],
  [Rainstorm], [1,294], [1,294], [1,294], [5,270], [5,270],
  [Severe], [413], [413], [413], [1,449], [1,449],
  table.cell(colspan: 6, fill: luma(235))[*Events detected (%) with 95% CI*],
  [All events], [44 #text(size: 6.5pt, fill: luma(90))[(39, 50)]], [71 #text(size: 6.5pt, fill: luma(90))[(67, 74)]], [73 #text(size: 6.5pt, fill: luma(90))[(70, 76)]], [66 #text(size: 6.5pt, fill: luma(90))[(63, 69)]], [72 #text(size: 6.5pt, fill: luma(90))[(70, 75)]],
  [Light], [35 #text(size: 6.5pt, fill: luma(90))[(30, 41)]], [63 #text(size: 6.5pt, fill: luma(90))[(59, 67)]], [67 #text(size: 6.5pt, fill: luma(90))[(63, 71)]], [54 #text(size: 6.5pt, fill: luma(90))[(51, 57)]], [62 #text(size: 6.5pt, fill: luma(90))[(59, 65)]],
  [Moderate], [58 #text(size: 6.5pt, fill: luma(90))[(50, 66)]], [87 #text(size: 6.5pt, fill: luma(90))[(84, 90)]], [87 #text(size: 6.5pt, fill: luma(90))[(83, 91)]], [88 #text(size: 6.5pt, fill: luma(90))[(86, 90)]], [90 #text(size: 6.5pt, fill: luma(90))[(88, 92)]],
  [Heavy], [72 #text(size: 6.5pt, fill: luma(90))[(64, 78)]], [89 #text(size: 6.5pt, fill: luma(90))[(85, 92)]], [90 #text(size: 6.5pt, fill: luma(90))[(86, 93)]], [94 #text(size: 6.5pt, fill: luma(90))[(92, 95)]], [95 #text(size: 6.5pt, fill: luma(90))[(93, 96)]],
  [Rainstorm], [74 #text(size: 6.5pt, fill: luma(90))[(67, 79)]], [91 #text(size: 6.5pt, fill: luma(90))[(88, 93)]], [91 #text(size: 6.5pt, fill: luma(90))[(87, 94)]], [93 #text(size: 6.5pt, fill: luma(90))[(92, 95)]], [94 #text(size: 6.5pt, fill: luma(90))[(93, 96)]],
  [Severe], [82 #text(size: 6.5pt, fill: luma(90))[(73, 89)]], [96 #text(size: 6.5pt, fill: luma(90))[(93, 99)]], [93 #text(size: 6.5pt, fill: luma(90))[(87, 96)]], [96 #text(size: 6.5pt, fill: luma(90))[(95, 98)]], [96 #text(size: 6.5pt, fill: luma(90))[(94, 98)]],
  table.cell(colspan: 6, fill: luma(235))[*Peak hour within ± 1 h (%)*],
  [All events], [49], [61], [61], [55], [55],
  [Light], [44], [58], [58], [54], [54],
  [Moderate], [55], [67], [66], [57], [57],
  [Heavy], [52], [64], [65], [54], [55],
  [Rainstorm], [59], [68], [70], [59], [58],
  [Severe], [70], [74], [80], [71], [72],
  table.cell(colspan: 6, fill: luma(235))[*Median rain-centre timing error (h, + = satellite late)*],
  [All events], [−0.52], [−0.21], [−0.26], [−0.33], [−0.38],
  [Light], [−0.82], [−0.39], [−0.44], [−0.52], [−0.53],
  [Moderate], [−0.62], [−0.20], [−0.23], [−0.38], [−0.35],
  [Heavy], [−0.21], [0.12], [0.03], [−0.09], [−0.17],
  [Rainstorm], [0.31], [0.31], [0.26], [0.12], [0.04],
  [Severe], [0.51], [0.41], [0.37], [0.44], [0.34],
  table.cell(colspan: 6, fill: luma(235))[*Best lag: late minus early share (percentage points)*],
  [All events], [−19], [−11], [−13], [−15], [−15],
  [Light], [−29], [−17], [−20], [−22], [−22],
  [Moderate], [−23], [−10], [−7], [−17], [−14],
  [Heavy], [−6], [6], [6], [−1], [−2],
  [Rainstorm], [14], [11], [12], [9], [7],
  [Severe], [34], [24], [23], [23], [18],
  table.cell(colspan: 6, fill: luma(235))[*Median satellite / gauge event duration*],
  [All events], [1.45], [1.75], [1.75], [1.25], [1.29],
  [Light], [2.00], [2.00], [2.00], [1.67], [1.67],
  [Moderate], [1.00], [1.33], [1.33], [1.12], [1.15],
  [Heavy], [1.00], [1.33], [1.24], [1.11], [1.09],
  [Rainstorm], [1.29], [1.50], [1.33], [1.16], [1.13],
  [Severe], [1.50], [1.54], [1.50], [1.33], [1.30],
  table.cell(colspan: 6, fill: luma(235))[*Pooled event volume bias (%) with 95% CI*],
  [All events], [6 #text(size: 6.5pt, fill: luma(90))[(−22, 41)]], [7 #text(size: 6.5pt, fill: luma(90))[(−4, 19)]], [32 #text(size: 6.5pt, fill: luma(90))[(10, 57)]], [−11 #text(size: 6.5pt, fill: luma(90))[(−16, −6)]], [13 #text(size: 6.5pt, fill: luma(90))[(1, 25)]],
  [Light], [133 #text(size: 6.5pt, fill: luma(90))[(54, 230)]], [125 #text(size: 6.5pt, fill: luma(90))[(94, 158)]], [198 #text(size: 6.5pt, fill: luma(90))[(141, 254)]], [36 #text(size: 6.5pt, fill: luma(90))[(21, 52)]], [83 #text(size: 6.5pt, fill: luma(90))[(58, 108)]],
  [Moderate], [−1 #text(size: 6.5pt, fill: luma(90))[(−32, 38)]], [24 #text(size: 6.5pt, fill: luma(90))[(7, 41)]], [47 #text(size: 6.5pt, fill: luma(90))[(26, 70)]], [−2 #text(size: 6.5pt, fill: luma(90))[(−10, 8)]], [31 #text(size: 6.5pt, fill: luma(90))[(14, 49)]],
  [Heavy], [−5 #text(size: 6.5pt, fill: luma(90))[(−40, 38)]], [−8 #text(size: 6.5pt, fill: luma(90))[(−20, 2)]], [19 #text(size: 6.5pt, fill: luma(90))[(−19, 64)]], [−8 #text(size: 6.5pt, fill: luma(90))[(−15, −2)]], [14 #text(size: 6.5pt, fill: luma(90))[(−4, 33)]],
  [Rainstorm], [−41 #text(size: 6.5pt, fill: luma(90))[(−62, −17)]], [−39 #text(size: 6.5pt, fill: luma(90))[(−48, −28)]], [−39 #text(size: 6.5pt, fill: luma(90))[(−54, −21)]], [−27 #text(size: 6.5pt, fill: luma(90))[(−34, −20)]], [−13 #text(size: 6.5pt, fill: luma(90))[(−30, 2)]],
  [Severe], [−54 #text(size: 6.5pt, fill: luma(90))[(−72, −35)]], [−62 #text(size: 6.5pt, fill: luma(90))[(−70, −54)]], [−62 #text(size: 6.5pt, fill: luma(90))[(−73, −51)]], [−56 #text(size: 6.5pt, fill: luma(90))[(−63, −49)]], [−51 #text(size: 6.5pt, fill: luma(90))[(−63, −38)]],
)

#let diurnal_rows = (
  table.cell(colspan: 6, fill: luma(235))[*Correlation of the diurnal amount curves*],
  [All seasons], [0.93], [0.91], [0.29], [0.90], [0.42],
  [MAM], [0.47], [0.69], [0.63], [0.88], [0.56],
  [JJA], [0.90], [0.97], [0.29], [0.97], [0.30],
  [SON], [0.57], [0.71], [0.34], [0.81], [0.68],
  [DJF], [0.06], [0.71], [0.81], [−0.17], [0.29],
  table.cell(colspan: 6, fill: luma(235))[*First-harmonic phase difference (h, + = satellite later)*],
  [All seasons], [−0.9], [−1.3], [−4.2], [−1.1], [−3.6],
  [MAM], [−3.4], [−2.4], [−2.6], [−1.0], [−2.8],
  [JJA], [1.4], [−0.3], [−4.2], [−0.5], [−3.9],
  [SON], [0.3], [0.1], [−1.0], [0.0], [−0.5],
  [DJF], [5.3], [2.1], [1.5], [5.4], [2.0],
  table.cell(colspan: 6, fill: luma(235))[*Gauge peak hour (BJT, hour-ending), over the same cells*],
  [All seasons], [18], [18], [18], [18], [18],
  [MAM], [15], [15], [15], [16], [16],
  [JJA], [18], [18], [18], [18], [18],
  [SON], [7], [7], [7], [11], [11],
  [DJF], [19], [19], [19], [13], [13],
  table.cell(colspan: 6, fill: luma(235))[*Satellite peak hour (BJT, hour-ending)*],
  [All seasons], [19], [17], [11], [17], [11],
  [MAM], [11], [15], [11], [16], [11],
  [JJA], [20], [17], [11], [17], [11],
  [SON], [10], [16], [11], [7], [11],
  [DJF], [22], [18], [19], [18], [22],
  table.cell(colspan: 6, fill: luma(235))[*Relative amplitude, satellite / gauge*],
  [All seasons], [1.52], [1.28], [1.77], [1.29], [1.94],
  [MAM], [1.19], [2.20], [1.45], [1.99], [2.31],
  [JJA], [1.72], [1.01], [1.19], [1.04], [1.16],
  [SON], [5.44], [3.58], [6.13], [2.30], [3.86],
  [DJF], [2.04], [2.01], [1.26], [0.43], [0.26],
)

#let series_rows = (
  [Median station r, hours either side is wet], [0.10], [0.32], [0.28], [0.34], [0.29],
  [Median station r, all hours (reference)], [0.23], [0.43], [0.37], [0.45], [0.39],
  [Regional-mean r, 1 h], [0.43], [0.84], [0.66], [0.83], [0.66],
  [Regional-mean r, 3 h], [0.43], [0.87], [0.68], [0.85], [0.70],
  [Regional-mean r, 24 h], [–], [–], [–], [0.91], [0.80],
  [Regional-mean KGE, 1 h], [−0.17], [0.79], [−0.04], [0.80], [0.01],
  [Regional-mean KGE, 24 h], [–], [–], [–], [0.89], [0.34],
  [Regional-mean relative bias, 1 h (%)], [28], [9], [30], [6], [33],
  [Regional-mean best lag, 1 h (h)], [0], [0], [0], [0], [0],
)

#let season_rows = (
  table.cell(colspan: 6, fill: luma(235))[*Relative bias, gauge-wet hours (%)*],
  [All seasons], [−48], [−30], [−11], [−32], [−11],
  [MAM], [−74], [−32], [−32], [−28], [−5],
  [JJA], [−62], [−31], [−14], [−31], [−16],
  [SON], [11], [−27], [18], [−30], [5],
  [DJF], [−94], [−65], [−55], [−82], [−80],
  table.cell(colspan: 6, fill: luma(235))[*Gauge-wet hours detected as rain (%)*],
  [All seasons], [34], [65], [66], [59], [60],
  [MAM], [23], [62], [62], [60], [62],
  [JJA], [38], [71], [68], [70], [68],
  [SON], [41], [65], [67], [64], [66],
  [DJF], [4], [24], [39], [14], [18],
  table.cell(colspan: 6, fill: luma(235))[*Events detected (%)*],
  [All seasons], [44], [71], [73], [66], [72],
  [MAM], [46], [77], [79], [74], [79],
  [JJA], [47], [75], [74], [79], [78],
  [SON], [49], [65], [73], [69], [79],
  [DJF], [1], [30], [51], [30], [42],
  table.cell(colspan: 6, fill: luma(235))[*Event peak within ± 1 h (%)*],
  [All seasons], [49], [61], [61], [55], [55],
  [MAM], [45], [64], [64], [58], [57],
  [JJA], [54], [62], [64], [58], [59],
  [SON], [41], [60], [57], [53], [53],
  [DJF], [14], [44], [41], [43], [44],
  table.cell(colspan: 6, fill: luma(235))[*Pooled event volume bias (%)*],
  [All seasons], [6], [7], [32], [−11], [13],
  [MAM], [−28], [16], [28], [−7], [23],
  [JJA], [−25], [2], [13], [−7], [5],
  [SON], [203], [32], [126], [−12], [27],
  [DJF], [−99], [−69], [−21], [−64], [−59],
  table.cell(colspan: 6, fill: luma(235))[*Median station r, hours either side is wet*],
  [All seasons], [0.10], [0.32], [0.28], [0.34], [0.29],
  [MAM], [0.02], [0.36], [0.32], [0.40], [0.33],
  [JJA], [0.09], [0.30], [0.26], [0.30], [0.26],
  [SON], [0.15], [0.42], [0.31], [0.40], [0.35],
  [DJF], [−0.08], [−0.10], [−0.03], [−0.20], [−0.17],
)

#let region_rows = (
  table.cell(colspan: 6, fill: luma(235))[*Relative bias, gauge-wet hours (%) with 95% CI*],
  [Relief < 500 m], [−49 #text(size: 6.5pt, fill: luma(90))[(−64, −31)]], [−29 #text(size: 6.5pt, fill: luma(90))[(−36, −22)]], [−12 #text(size: 6.5pt, fill: luma(90))[(−30, 6)]], [−32 #text(size: 6.5pt, fill: luma(90))[(−38, −27)]], [−11 #text(size: 6.5pt, fill: luma(90))[(−25, 2)]],
  [Relief 500–1000 m], [−47 #text(size: 6.5pt, fill: luma(90))[(−60, −31)]], [−31 #text(size: 6.5pt, fill: luma(90))[(−37, −25)]], [−8 #text(size: 6.5pt, fill: luma(90))[(−23, 9)]], [−32 #text(size: 6.5pt, fill: luma(90))[(−37, −27)]], [−8 #text(size: 6.5pt, fill: luma(90))[(−20, 2)]],
  [#text(fill: luma(110))[Relief ≥ 1000 m#super[a]]], [−54 #text(size: 6.5pt, fill: luma(90))[(−66, −41)]], [−38 #text(size: 6.5pt, fill: luma(90))[(−45, −31)]], [−18 #text(size: 6.5pt, fill: luma(90))[(−33, −4)]], [−37 #text(size: 6.5pt, fill: luma(90))[(−43, −31)]], [−18 #text(size: 6.5pt, fill: luma(90))[(−28, −7)]],
  table.cell(colspan: 6, fill: luma(235))[*Relative bias, heavy hours 4–8 mm/h (%) with 95% CI*],
  [Relief < 500 m], [−58 #text(size: 6.5pt, fill: luma(90))[(−71, −39)]], [−44 #text(size: 6.5pt, fill: luma(90))[(−52, −36)]], [−20 #text(size: 6.5pt, fill: luma(90))[(−47, 9)]], [−42 #text(size: 6.5pt, fill: luma(90))[(−48, −36)]], [−18 #text(size: 6.5pt, fill: luma(90))[(−38, 2)]],
  [Relief 500–1000 m], [−55 #text(size: 6.5pt, fill: luma(90))[(−67, −42)]], [−39 #text(size: 6.5pt, fill: luma(90))[(−46, −30)]], [−16 #text(size: 6.5pt, fill: luma(90))[(−36, 3)]], [−33 #text(size: 6.5pt, fill: luma(90))[(−41, −26)]], [−16 #text(size: 6.5pt, fill: luma(90))[(−30, −1)]],
  [#text(fill: luma(110))[Relief ≥ 1000 m#super[a]]], [−49 #text(size: 6.5pt, fill: luma(90))[(−65, −29)]], [−39 #text(size: 6.5pt, fill: luma(90))[(−48, −28)]], [−10 #text(size: 6.5pt, fill: luma(90))[(−33, 17)]], [−33 #text(size: 6.5pt, fill: luma(90))[(−44, −21)]], [−14 #text(size: 6.5pt, fill: luma(90))[(−31, 6)]],
  table.cell(colspan: 6, fill: luma(235))[*Events detected (%) with 95% CI*],
  [Relief < 500 m], [46 #text(size: 6.5pt, fill: luma(90))[(40, 51)]], [72 #text(size: 6.5pt, fill: luma(90))[(68, 76)]], [74 #text(size: 6.5pt, fill: luma(90))[(70, 78)]], [68 #text(size: 6.5pt, fill: luma(90))[(65, 71)]], [73 #text(size: 6.5pt, fill: luma(90))[(70, 76)]],
  [Relief 500–1000 m], [43 #text(size: 6.5pt, fill: luma(90))[(38, 48)]], [69 #text(size: 6.5pt, fill: luma(90))[(65, 72)]], [73 #text(size: 6.5pt, fill: luma(90))[(69, 76)]], [65 #text(size: 6.5pt, fill: luma(90))[(62, 68)]], [72 #text(size: 6.5pt, fill: luma(90))[(69, 74)]],
  [#text(fill: luma(110))[Relief ≥ 1000 m#super[a]]], [42 #text(size: 6.5pt, fill: luma(90))[(36, 48)]], [67 #text(size: 6.5pt, fill: luma(90))[(63, 71)]], [71 #text(size: 6.5pt, fill: luma(90))[(66, 75)]], [63 #text(size: 6.5pt, fill: luma(90))[(60, 66)]], [71 #text(size: 6.5pt, fill: luma(90))[(68, 74)]],
  table.cell(colspan: 6, fill: luma(235))[*Event peak within ± 1 h (%)*],
  [Relief < 500 m], [50], [61], [61], [55], [56],
  [Relief 500–1000 m], [49], [63], [62], [56], [55],
  [#text(fill: luma(110))[Relief ≥ 1000 m#super[a]]], [44], [60], [58], [53], [54],
  table.cell(colspan: 6, fill: luma(235))[*Median station r, hours either side is wet*],
  [Relief < 500 m], [0.10], [0.32], [0.27], [0.34], [0.29],
  [Relief 500–1000 m], [0.09], [0.31], [0.29], [0.34], [0.30],
  [#text(fill: luma(110))[Relief ≥ 1000 m#super[a]]], [0.10], [0.36], [0.33], [0.39], [0.34],
  table.cell(colspan: 6, fill: luma(235))[*Diurnal phase difference (h, + = satellite later)*],
  [Relief < 500 m], [−1.3], [−1.6], [−4.4], [−1.2], [−4.1],
  [Relief 500–1000 m], [−0.7], [−1.3], [−4.1], [−1.2], [−3.3],
  [#text(fill: luma(110))[Relief ≥ 1000 m#super[a]]], [−1.5], [−1.2], [−3.3], [−1.4], [−2.6],
)

#let note(body) = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, fill: luma(246),
  stroke: (left: 2.5pt + rgb("#2a78d6")), body)

// ---------------------------------------------------------------------------------------------
// Title

#align(left)[
  #set text(hyphenate: false)
  #text(font: ("Segoe UI", "Arial"), size: 9pt, fill: luma(100), tracking: 0.04em)[INTERIM RESULTS · 13 SEPTEMBER 2026]
  #v(0.2em)
  #text(font: ("Segoe UI", "Arial"), size: 19pt, weight: "semibold")[Can satellite precipitation products capture the temporal evolution of rainfall?]
  #v(0.1em)
  #text(size: 12pt, fill: luma(60))[FY4B, GPM IMERG and GSMaP against 237 hourly gauges in the Shiyan study area, 2022–2024 — by rain intensity, season and landform region]
]
#v(0.6em)
#line(length: 100%, stroke: 0.5pt + luma(180))

= Key findings

- *GPM follows the hour-to-hour evolution of rainfall best, GSMaP second, FY4B last.* Averaged over the gauge network, the hourly correlation is 0.83 (GPM), 0.66 (GSMaP) and 0.43 (FY4B, on the hours it covers). At single gauges the median correlation over rainy hours is only 0.34, 0.29 and 0.10: point-scale hourly timing is weak for every product.
- *Rain events are usually detected, but timing is coarse.* GPM detects 66% of gauge rain events and GSMaP 72%; the event peak falls within ± 1 h of the gauge peak in 55% of events for both. Satellite rain is centred slightly *early* in light events (about −0.5 h) and *late* in the most intense events (+0.3 to +0.4 h).
- *Underestimation grows steadily with intensity.* GPM is 11% too low on light hours and 85% too low on hours ≥ 20 mm/h (GSMaP: +16% and −80%). The most intense hours are nearly always detected as rain (85–89%) but almost never in the right class (1–6%).
- *Winter (DJF) is where the products fail.* Hourly bias drops to about −80%, only 14–18% of rainy gauge hours are detected, and FY4B detects 1% of winter events. Part of this is a gauge artefact: the gauges' winter diurnal cycle peaks at midday, which is the signature of delayed snowmelt tips.
- *Landform makes no measurable difference.* No gauge sits on a plain; between the two relief regions with enough gauges (< 500 m and 500–1000 m) every metric agrees within its confidence interval.
- *GSMaP carries a spurious late-morning maximum* (09–12 BJT, three times the gauge amount at 11:00 in summer) that shifts its diurnal phase by about −4 h. It is not a clock offset: the best lag of its regional series is 0 h.

= Question and scope

The question has four parts: whether the products reproduce the temporal evolution of rainfall at all, and whether that ability changes with rain intensity, with season and with geographical region. The earlier heavy-rain event study addressed the *spatial* pattern of daily totals on selected days; this study addresses *time*, over every hour of the record. Hours when the gauge reports no rain are excluded from the intensity analysis, as requested, and reported separately as false alarms (@sec-false-alarms), because a product that rains all the time would otherwise look perfect.

= Data

- *Gauges.* 237 stations in the study area (109.58–111.47° E, 31.25–33.20° N), hourly, labelled hour-ending Beijing time (BJT). They are *tipping buckets with 0.5 mm resolution*: 99.9% of rainy values in 2022–2024 are multiples of 0.5 mm.
- *Satellite products*, sampled at the gauge pixels and aligned to the gauge clock: FY4B precipitation (strict hourly aggregation, navigation-corrected); GPM IMERG V07 (Google Earth Engine collection `NASA/GPM_L3/IMERG_V07`); GSMaP V8 operational `hourlyPrecipRate`, which carries no gauge correction.
- *Period.* The 08–08 BJT days from 1 January 2022 to 31 December 2024: 26,304 hours, contiguous.
- *Two evaluation samples* (@tab-coverage). FY4B never provides the 23:00 label and misses whole months (2022-01 to 05, 2022-10 to 2023-03, 2024-02). A fair three-product comparison is therefore restricted to the hours FY4B covers; GPM and GSMaP are also scored over the full record, which is the basis for the seasonal analysis. Within each sample every product is scored on identical station-hours.

#block(breakable: false)[
#figure(
  table(
    columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
    align: (left, right, right, right, right, right, right, right, right),
    stroke: none,
    inset: (x: 4pt, y: 2.6pt),
    table.vline(x: 5, stroke: 0.4pt + luma(170)),
    table.header(
      rule,
      table.cell(rowspan: 2, align: left + bottom)[*Season*],
      table.cell(colspan: 4, align: center)[*Hours FY4B covers* (three products)],
      table.cell(colspan: 4, align: center)[*Full record* (GPM, GSMaP)],
      [*Hours*], [*Days*], [*Wet*#super[b]], [*Wet %*],
      [*Hours*], [*Days*], [*Wet*#super[b]], [*Wet %*],
      thin,
    ),
    ..coverage_rows,
    rule,
  ),
  caption: [Evaluation samples. #super[b] Wet station-hours: gauge rain ≥ 0.1 mm; the wet share is relative to all station-hours scored. Seasons follow the month of the 08–08 BJT day.],
) <tab-coverage>
]

= Methods

== Intensity classes, seasons and uncertainty

#block(breakable: false)[
#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, center, left),
    stroke: none,
    inset: (x: 6pt, y: 2.6pt),
    table.header(rule, [*Class*], [*Hourly precipitation*], [*Treatment*], thin),
    [No rain / trace], [< 0.1 mm/h], [excluded from the class metrics; reported as false alarms],
    [Short-duration light rain], [\[0.1, 2) mm/h], [for gauges this means 0.5, 1.0 or 1.5 mm],
    [Short-duration moderate rain], [\[2, 4) mm/h], [],
    [Short-duration heavy rain], [\[4, 8) mm/h], [],
    [Short-duration rainstorm], [\[8, 20) mm/h], [],
    [Severe rainstorm and above], [≥ 20 mm/h], [],
    rule,
  ),
  caption: [Hourly intensity classes (lower bound inclusive). Hours are classified by the gauge value; the satellite value is classified with the same bounds to test class agreement.],
) <tab-classes>
]

Seasons are meteorological: MAM, JJA, SON and DJF. Uncertainty is a *day-block bootstrap*: whole 08–08 days are resampled with replacement (1000 replicates, fixed seed), which keeps the hours and stations of one weather system together. Intervals are given for additive statistics (bias, detection rates, class agreement, pooled event volume bias); medians are reported without intervals. A stratum is flagged as low-sample below 100 station-hours, 10 days or 30 events.

== Landform regions <sec-landform-method>

The requested scheme is a plain–hills–mountains regionalization. No geomorphological boundary data exist in the repository, so a reproducible rule was derived from the 30 m Copernicus GLO-30 DEM:

+ *Classifier: local relief amplitude* — maximum minus minimum elevation in a square window centred on the gauge — with the breaks of the Chinese basic landform scheme: plain < 30 m, hills 30–200 m, and mountains split at 500 and 1000 m into small-, medium- and large-relief mountains.
+ *Window fixed before classification* by the mean change-point method on the area-wide mean relief of non-overlapping windows from 0.3 to 9 km: the change point of ln(relief / area) falls at *2.7 km* (@fig-landform, right).
+ *Slope and elevation do not reassign gauges.* Mean window slope rises monotonically with relief class (@tab-landform), which confirms the classes are physically coherent.
+ A region needs at least *20 gauges* to be interpreted.

#note[
  *Limitation.* Under these breaks the network contains no plain and only 13 hills, so a true plain–hills–mountains comparison is not possible. The evaluation regions are the relief groups < 500 m (129 gauges) and 500–1000 m (97); ≥ 1000 m (11 gauges) is reported but not interpreted. Absolute labels depend on the window — 160 gauges change region at half or double the window — but the ranking of gauges by ruggedness does not (Spearman 0.91–0.92), so a comparison between less and more rugged gauges is robust. Ruggedness also increases southward, so regional differences cannot be attributed to terrain alone. The class breaks match those of a published regionalization (Sci. Rep., 2020, doi:10.1038/s41598-020-66993-9); the change-point sequence ln(relief / area) follows the usual description of the method and still has to be checked against its original source before publication.
]

== Metrics

- *Hourly, by gauge class* (gauge-wet hours): relative bias of the amount, detection as rain (satellite ≥ 0.1 mm), share of hours in the same class and in a lower class, and the full gauge-class × satellite-class confusion matrix.
- *False alarms* (gauge-dry hours): how often the satellite reports ≥ 0.1 mm and ≥ 0.5 mm, the false-alarm ratio, and the share of satellite rain volume falling on those hours.
- *Station rain events.* At each gauge an event is a run of wet hours ending after ≥ 3 dry hours (3 h because a 0.5 mm bucket tips only every few hours in drizzle; a 1 h gap was run as a sensitivity). The event class is its peak hourly intensity. Each product is scored over the event ± 3 h, cut back so the window never reaches another event, and only where gauge and products are complete over that window. Scores: detection, peak and rain-centre timing error, onset and end error, best lag within ± 3 h, duration ratio, volume bias and peak ratio.
- *Diurnal cycle*: mean amount and wet-hour frequency per hour of day; correlation of the curves and first-harmonic phase difference.
- *Hourly correlation at each gauge*, over the hours where the gauge *or* the satellite is wet, so that the shared dry hours do not inflate it (the all-hours value is kept as a reference).
- *Regional-mean series*: hourly means over each region's gauges (kept when ≥ 90% report), accumulated to 3 h and 24 h; correlation, KGE, bias and best lag.

= Results

== Landform regions

#figure(
  image("figures/fig1_landform_regions.png", width: 100%),
  caption: [Local relief in 2.7 km windows with the gauges by region (left); area-wide mean relief and the change-point gain against window size (right).],
) <fig-landform>

#block(breakable: false)[
#figure(
  table(
    columns: (1fr, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    stroke: none,
    inset: (x: 6pt, y: 2.6pt),
    table.header(rule, [*Class*], [*Gauges*], [*Mean relief (m)*], [*Mean window slope (°)*], [*Mean elevation (m)*], thin),
    ..landform_rows,
    rule,
  ),
  caption: [Gauge counts by landform class at the 2.7 km window. #super[a] Below the 20-gauge minimum: reported, not interpreted.],
) <tab-landform>
]

== Hourly performance by intensity class

Every product underestimates more as the rain gets heavier (@tab-intensity, @fig-intensity). On light hours GPM is within about 10% of the gauges and GSMaP is too wet; from 8 mm/h upward all products lose roughly 60–86% of the gauge amount. Detection behaves differently: the most intense hours are the ones most reliably *detected* as rain (85–89% for GPM and GSMaP), yet only 1–6% are placed in the ≥ 20 mm/h class — the satellites see the storm but smear its intensity. The confusion matrices (@fig-confusion) show where the hours go: most gauge rainstorm and severe hours land in the satellite's light or moderate classes. FY4B misses two-thirds of all rainy gauge hours and is the weakest at every intensity.

#figure(
  split-table("Gauge class (mm/h)", intensity_rows),
  caption: [Hourly metrics by gauge intensity class, all seasons and gauges. Grey figures in brackets: 95% day-block bootstrap interval.],
) <tab-intensity>

#figure(
  image("figures/fig2_intensity_classes_gpm_gsmap_full.png", width: 100%),
  caption: [Relative bias, detection, class agreement and underestimation by gauge class and season, GPM and GSMaP over the full record. The FY4B-hours version is in @fig-intensity-fy4b.],
) <fig-intensity>

#figure(
  image("figures/fig3_intensity_confusion.png", width: 94%),
  caption: [Where each gauge class lands in the satellite classes (row percentages). Outlined cells: the same class.],
) <fig-confusion>

== False alarms <sec-false-alarms>

On 5.6–7.2% of gauge-dry hours the satellites report rain, and about 60% of all satellite-wet hours are gauge-dry (70% for FY4B). Much of this is light: raising the satellite threshold to the gauge resolution of 0.5 mm cuts the rate by more than half for GPM and GSMaP, and by about a third for FY4B (@tab-false-alarms, @fig-false-alarms). The share of satellite rain *volume* falling on gauge-dry hours is about a third for GPM and GSMaP but 60% for FY4B, so FY4B's rain is largely displaced in time or space. Some of the 0.1–0.5 mm false alarms may be real drizzle that the bucket had not yet tipped for.

#figure(
  split-table("Season", false_alarm_rows, first-width: 18%),
  caption: [False alarms on gauge-dry hours (gauge < 0.1 mm), all gauges.],
) <tab-false-alarms>

#figure(
  image("figures/fig4_false_alarms.png", width: 100%),
  caption: [Satellite rain on gauge-dry hours by season, split by the satellite's hourly amount.],
) <fig-false-alarms>

== Rain events: detection, timing and volume

Across 86,508 complete gauge events (full record), detection rises from 54–62% for light events to 96% for severe ones, and the share of peaks within ± 1 h rises from 54% to 71–72% (@tab-events, @fig-events). At hourly resolution the median peak-hour error is 0 h almost everywhere, so the sub-hour structure is read from the rain-weighted centre of each event. It shows a consistent drift with intensity for all three products: satellite rain is centred about 0.4–0.8 h *early* in light events and 0.3–0.5 h *late* in severe ones, and the best-lag balance moves from early to late in step. Satellite events last longer than gauge events (median ratio 1.1–2.0), light events carry far too much rain (+36% GPM, +83% GSMaP) and severe events far too little (−56%, −51%).

Two caveats bear on light events: the 0.5 mm bucket delays the recorded onset of drizzle, which contributes to the apparent early satellite onset, and FY4B's gaps leave the three-product sample only 28% of all gauge events, biased toward short events that avoid 23:00.

#figure(
  split-table("Event class", event_rows, first-width: 20%),
  caption: [Event scores, all seasons and gauges. Timing statistics are over detected events; brackets: 95% day-block bootstrap interval.],
) <tab-events>

#figure(
  image("figures/fig5_event_timing.png", width: 100%),
  caption: [Event detection, peak timing, rain-centre timing error, best-lag balance and volume bias by event class.],
) <fig-events>

== Diurnal cycle

In summer the gauges peak at 18:00 BJT. GPM reproduces the shape almost exactly (r = 0.97, phase −0.5 h), FY4B peaks late (+1.4 h) and GSMaP peaks at 11:00 instead (@tab-diurnal, @fig-diurnal). The GSMaP maximum at 09–12 BJT appears in spring, summer and autumn and in both relief regions; its amount at 11:00 in summer is 3.0 times the gauge's while its wet-hour frequency rises much less, so it is a burst of excessive amounts rather than a clock shift (its regional series has best lag 0 h, @tab-series). In winter the gauge cycle itself peaks at 13:00, which fits snowmelt tips rather than rainfall timing, so the large winter phase differences are not interpretable as satellite errors.

#figure(
  split-table("Season", diurnal_rows, first-width: 18%),
  caption: [Diurnal cycle of the mean hourly amount, all gauges. In the FY4B-hours sample, SON covers only September 2022, September 2023 and September–November 2024, and DJF only January and December 2024.],
) <tab-diurnal>

#figure(
  image("figures/fig6_diurnal_season_gpm_gsmap_full.png", width: 84%),
  caption: [Diurnal cycles of mean amount and wet-hour frequency by season, GPM and GSMaP over the full record. The three-product version is @fig-diurnal-fy4b.],
) <fig-diurnal>

== Hourly agreement at gauges and in regional series

Point-scale hourly agreement is low for all products (@tab-series, @fig-stations): the median gauge correlation over rainy hours is 0.34 for GPM and 0.29 for GSMaP. Averaging over the network removes much of the point noise, raising GPM to 0.83 at 1 h and 0.91 at 24 h with KGE up to 0.89. GSMaP's regional series correlate reasonably (0.66–0.80) but its KGE stays low because of a +33% volume bias. No product shows a systematic hourly offset: the best lag of every regional series is 0 h.

#block(breakable: false)[
#figure(
  split-table("Statistic (all gauges, all seasons)", series_rows, first-width: 34%),
  caption: [Hourly agreement at gauges and of network-mean series. FY4B never completes a 24-hour day, so its sample has no 24 h values.],
) <tab-series>
]

#figure(
  image("figures/fig7_station_correlation_maps.png", width: 88%),
  caption: [Hourly correlation at each gauge over hours where gauge or satellite is wet; marker shape shows the relief region.],
) <fig-stations>

#figure(
  image("figures/fig9_regional_series.png", width: 92%),
  caption: [Correlation, KGE and relative bias of regional-mean series at 1, 3 and 24 h.],
) <fig-series>

== Seasonal differences

Spring, summer and autumn are broadly similar for GPM and GSMaP (@tab-seasons, @fig-scorecard): hourly bias −28% to −31% for GPM, event detection 69–79%, peaks within ± 1 h in 53–59% of events. Winter stands apart on every metric — hourly bias −82% (GPM) and −80% (GSMaP), 14–18% of rainy gauge hours detected, 30–42% of events detected, negative station correlations — and FY4B detects almost nothing in winter. Because the winter gauge record is partly snowmelt, winter numbers mix satellite failure (snowfall and shallow cold-season rain are hard to retrieve) with gauge timing error. FY4B also varies strongly outside winter: −74% in spring but an overestimate in autumn, where its event volume is three times the gauge's.

#figure(
  split-table("Season", season_rows, first-width: 18%),
  caption: [Key metrics by season, all gauges.],
) <tab-seasons>

#figure(
  image("figures/fig8_scorecard_gpm_gsmap_full.png", width: 80%),
  caption: [Scorecard by season and landform region, GPM and GSMaP over the full record. The three-product version is @fig-scorecard-fy4b.],
) <fig-scorecard>

== Regional differences

The two interpretable relief regions are statistically indistinguishable (@tab-regions, @fig-diurnal-region). GPM's hourly bias is −32% in both, its event detection 68% against 65% and its median station correlation 0.34 in both; for heavy hours the difference (−42% against −33%) lies well inside overlapping intervals. The ≥ 1000 m group points the same way but has too few gauges to interpret. Within this network, therefore, the temporal skill of the products is set by rain intensity and season, not by local terrain ruggedness.

#figure(
  split-table("Region", region_rows, first-width: 22%),
  caption: [Key metrics by relief region, all seasons. #super[a] 11 gauges, below the 20-gauge minimum: not interpreted.],
) <tab-regions>

#figure(
  image("figures/fig6_diurnal_region_gpm_gsmap_full.png", width: 84%),
  caption: [Diurnal cycles by relief region, GPM and GSMaP over the full record.],
) <fig-diurnal-region>

= Answers to the research questions

+ *Do the products capture the temporal evolution of rainfall?* Partly, and at the network scale far better than at a gauge. GPM tracks the regional hourly series closely (r = 0.83, KGE 0.80) and reproduces the summer diurnal cycle; at single gauges hourly correlation is 0.3 at best, event peaks are within ± 1 h only about half the time, and satellite events are longer and smoother than gauge events. GSMaP detects events slightly more often but has a spurious late-morning maximum and a wet volume bias. FY4B is the weakest on every temporal measure.
+ *Does performance depend on intensity?* Strongly. Detection and peak timing improve with intensity, while the amount is progressively underestimated (−80% to −86% at ≥ 20 mm/h) and class agreement collapses. Light rain is overestimated in event volume and duration. Timing drifts from slightly early (light) to slightly late (intense).
+ *Does it depend on season?* Winter is by far the worst season for all products; spring, summer and autumn differ little for GPM and GSMaP. Winter results are partly contaminated by snowmelt timing in the gauge record.
+ *Does it depend on region?* Not detectably between the two relief regions the network supports. The network holds no plains and too few gauges in the most rugged terrain to test those extremes.

= Limitations and open issues

- *Gauge resolution.* 0.5 mm tipping buckets cannot resolve 0.1–0.5 mm/h and delay drizzle onsets; light-rain detection, false alarms and early satellite onsets are affected.
- *Winter gauge record.* The midday winter gauge peak suggests delayed snowmelt tips; winter timing metrics should not be read as satellite errors without a snow screen.
- *FY4B coverage.* The three-product sample covers 13,243 hours, 28% of gauge events, and winter only in January and December 2024; FY4B seasonal values rest on uneven months.
- *Landform.* No plains; region labels depend on the relief window (ranking is stable); terrain is confounded with latitude; the change-point sequence definition awaits verification against its original source, which could not be retrieved for this report.
- *GSMaP late-morning maximum.* Cause unknown. It should be checked against the raw Earth Engine export and another GSMaP product version before being attributed to the retrieval.
- *Point versus pixel.* All comparisons are gauge against the pixel containing it; part of the low station-scale agreement is representativeness error, not retrieval error.

= Next steps

+ Verify the mean change-point sequence definition against the original method papers on optimal relief windows.
+ Screen winter hours for snow — for example with the ERA5-Land 2 m temperature already extracted at the gauges — and rerun the DJF metrics.
+ Inspect the raw GSMaP export for the 09–12 BJT maximum and compare with GSMaP gauge-corrected or reanalysis versions.
+ Repeat the event analysis with a dry-gap of 1 h (already computed) and 6 h to confirm the timing drift is not a segmentation artefact.

#pagebreak()

#heading(numbering: none)[Appendix A — Three-product figures]

#figure(
  image("figures/fig2_intensity_classes_all_products.png", width: 100%),
  caption: [Hourly performance by gauge class and season on the hours FY4B covers.],
) <fig-intensity-fy4b>

#figure(
  image("figures/fig6_diurnal_season_all_products.png", width: 84%),
  caption: [Diurnal cycles by season on the hours FY4B covers (no 23:00 values).],
) <fig-diurnal-fy4b>

#figure(
  image("figures/fig8_scorecard_all_products.png", width: 80%),
  caption: [Scorecard by season and landform region on the hours FY4B covers.],
) <fig-scorecard-fy4b>

#heading(numbering: none)[Appendix B — Reproduction]

All results come from commit `60012c4` on branch `satellite-temporal-evaluation`.

```sh
julia --project=. scripts/prepare_station_landform.jl          # landform regions (GDAL + GMT)
julia --project=. scripts/run_satellite_temporal_evaluation.jl # all CSVs, about 1.5 min
py -3.13 scripts/plot_satellite_temporal_evaluation.py         # figures
```

#table(
  columns: (58%, 42%),
  stroke: none,
  inset: (x: 5pt, y: 2.4pt),
  table.header(rule, [*File*], [*Content*], thin),
  [`data/processed/covariates/station_landform.csv`], [relief, class and region of every gauge],
  [`output/satellite_temporal_evaluation/landform/`], [window change point, class counts, window sensitivity, relief grid],
  [`sample_coverage.csv`], [hours, days and wet station-hours per sample, season and region],
  [`intensity_class_metrics.csv`], [hourly metrics by gauge class with bootstrap intervals],
  [`intensity_confusion.csv`], [gauge class × satellite class counts],
  [`false_alarm_metrics.csv`], [gauge-dry hours],
  [`station_events.csv` \ `event_scores.csv`], [every gauge event and its per-product scores],
  [`event_timing_summary.csv`], [event metrics by class, season and region (dry gaps 3 h and 1 h)],
  [`diurnal_cycle.csv` \ `diurnal_summary.csv`], [composite cycles and their phase statistics],
  [`station_hourly_correlation.csv` \ `station_correlation_summary.csv`], [gauge-level correlation],
  [`regional_series_metrics.csv`], [regional-mean series at 1, 3 and 24 h],
  [`run_settings.csv`], [every threshold, window and seed],
  rule,
)
