"""
What the benchmark's satellite inputs lost when the gridded products were reduced to gauge points.

The benchmark never opens a satellite grid: every product reaches it as a wide `time x station_id`
table, sampled nearest-pixel from the FY4B disk (`FY4BPreprocessing.extract_precipitation`) or as a
0.1 degree cell mean from Earth Engine. This module measures the cost of that reduction, and it
measures only - it fits nothing, changes no benchmark output, and is not part of any run.

Six questions, one group of functions each:

- `pixel_assignment_table` / `effective_cell_keys`: which grid cell each gauge actually reads, on
  the FY4B disk and on the 0.1 degree lat/lon grid. FY4B needs the second function because the
  archive holds segments from two sub-satellite longitudes, so one gauge has two pixels.
- `cell_membership_table` / `cell_multiplicity_summary`: how many gauges share a cell, and so how
  many receive a satellite value that carries no information about where inside the cell they are.
- `identical_series_groups` / `identical_series_table` / `cell_grouping_agreement`: the same
  question answered from the shipped tables alone, with no grid geometry assumed, and then used to
  test the geometry. Both assumptions needed it: the single-projection FY4B grouping and the
  0.1 degree grid the Earth Engine exports were assumed to sample are each contradicted by what
  the stored values can actually distinguish. Where the two disagree, the data is the answer.
- `gauge_pair_table` / `pair_disagreement_summary`: how much two gauges inside one cell disagree
  with each other. That is the floor on point-to-pixel error - a satellite cannot match both - and
  the ceiling on what any downscaling could recover. Pairs that do *not* share a cell, at the same
  separation, are the control.
- `pixel_footprint_km` / `footprint_table`: FY4B's 4 km is at nadir, and the study area is 23
  degrees of longitude and 32 degrees of latitude away from the sub-satellite point. This is the
  real ground footprint there.
- `domain_raster_values` / `terrain_representativeness_table`: how biased the gauge network is as a
  sample of the terrain it is used to validate over.

Conventions shared with the rest of the benchmark: matrices are `[station, time]` with `NaN` for
missing; a wet hour is `>= WET_MM`; hour labels are hour-ending Beijing time.
"""
module GridSupportDiagnostics

using DataFrames, Dates, Statistics
using NCDatasets
using MixedGWR: metric_continuous
using Main: FY4BPreprocessing
using Main.TraditionalInterpolation: haversine_distance_matrix

export WET_MM, COARSE_CELL_DEG, FY4B_NADIR_RESOLUTION_M, SEPARATION_BIN_EDGES
export fy4b_pixel_index, latlon_cell_index, latlon_cell_center, pixel_footprint_km
export pixel_assignment_table, effective_cell_keys
export cell_membership_table, cell_multiplicity_summary
export identical_series_groups, identical_series_table
export series_group_keys, cell_grouping_agreement
export pair_disagreement, gauge_pair_table, pair_disagreement_summary
export footprint_table, domain_raster_values, terrain_representativeness_table
export bilinear_neighbours, bilinear_combine, sample_precipitation_modes
export fy4b_hourly_modes, extraction_sensitivity_table

"""A wet hour, matching the benchmark's `rain_threshold` and `wet_threshold`."""
const WET_MM = 0.1

"""
The lat/lon cell size GPM and GSMaP were assumed to be sampled on: `scale: 11132` m in the Earth
Engine exports, which is 0.1 degree at the equator. An assumption - see `cell_grouping_agreement`.
"""
const COARSE_CELL_DEG = 0.1

"""FY4B's nominal resolution, which `FY4BPreprocessing.RESOLUTION` states is *at nadir*."""
const FY4B_NADIR_RESOLUTION_M = 4000.0

"""
Earth radius for converting a grid step to ground km.

`FY4BPreprocessing`'s value rather than `TraditionalInterpolation`'s 6378.388: the footprint is a
property of the geostationary projection, so it has to use the same figure that projection does.
The two differ by 0.25 km in 6378, far below anything reported here.
"""
const EARTH_RADIUS_KM = 6378.137

"""Gauge separations, in km, that `pair_disagreement_summary` bins pairs into."""
const SEPARATION_BIN_EDGES = [0.0, 5.0, 10.0, 20.0, 40.0, 60.0]

include("grid_support/cell_index.jl")
include("grid_support/footprint_geometry.jl")
include("grid_support/cell_membership.jl")
include("grid_support/identical_series.jl")
include("grid_support/pair_disagreement.jl")
include("grid_support/representativeness.jl")
include("grid_support/extraction_sensitivity.jl")

end # module
