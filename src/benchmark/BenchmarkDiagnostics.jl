"""
Post-hoc diagnostics for `run_interpolation_benchmark` results.

Everything here works off the artefacts an existing benchmark run already wrote
(`oof_*.csv`, `common_evaluation_mask.csv`, `split_common.csv`, `parameter_scan.csv`, ...),
so a diagnosis costs no benchmark re-run.

The scoring helper `_metrics` is a deliberate independent reimplementation of
`MixedGWR.metric_continuous` rather than a call into it: reproducing `metrics_pooled.csv`
from these matrices is the correctness check for the whole reading path, and sharing the
scoring code would make that check vacuous.

That is now the only reason. `metric_continuous` used to live in `MGERPipeline.jl`, a top-level
script a module could not import from at all; it has since moved into `MixedGWR`, so the
duplication here is a choice rather than a constraint. The two also differ at the edges -
`_metrics` returns a NaN row for an empty sample where `metric_continuous` asserts, and adds
`MSE`/`variance` - so they are not interchangeable as written.
"""
module BenchmarkDiagnostics

using CSV, DataFrames, Dates, Statistics

export read_wide_matrix, read_mask_matrix, read_fold_map, load_prediction_matrices
export null_baseline_matrices, satellite_rescale_matrix
export null_baseline_table, satellite_offset_table, mse_decomposition_table
export bandwidth_saturation_table, covariate_contribution_table
export rebuild_common_mask, dropout_table, mask_cost_table, run_comparison_table
export benchmark_diagnostics_outdir
export RAIN_CLASSES, MASK_METHODS

"""Rain-intensity strata, mirroring `InterpolationBenchmark.append_stratified_metrics!`."""
const RAIN_CLASSES = [
    ("no_rain", -Inf, 0.1), ("light", 0.1, 2.5),
    ("moderate", 2.5, 8.0), ("heavy", 8.0, Inf),
]

"""
Methods whose finiteness defines the shared evaluation mask, mirroring
`InterpolationBenchmark.MASK_METHODS`. Duplicated rather than imported because
`InterpolationBenchmark.jl` is a top-level script, not a module; `rebuild_common_mask` is
checked against the run's own `common_evaluation_mask.csv` so the copy cannot drift silently.
"""
const MASK_METHODS = [
    "raw", "idw", "adw", "tps", "gwr", "residual_gwr", "mixed_gwr", "mgwr",
]

# The diagnostics themselves are split into one file per question, in `diagnostics/`, keeping the
# `D<n>` labels the sections have always carried. They are `include`d into this module rather than
# being modules of their own, so every name below stays exactly where callers already reach it -
# including the private `_metrics` and `_rain_mask` that `test-benchmark-diagnostics.jl` calls
# qualified. Only `reading` and `scoring` are shared; the `D` files do not call each other.

include("diagnostics/reading.jl")
include("diagnostics/scoring.jl")
include("diagnostics/null_baselines.jl")
include("diagnostics/satellite_offset.jl")
include("diagnostics/mse_decomposition.jl")
include("diagnostics/mask_cost.jl")
include("diagnostics/run_comparison.jl")
include("diagnostics/bandwidth_saturation.jl")
include("diagnostics/covariate_contribution.jl")

end # module
