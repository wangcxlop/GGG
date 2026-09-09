"""
The covariate design vocabulary: which groups exist, which columns each is made of, and which
family a group belongs to.

`JointCovariateModels` (which fits the designs) and `JointVariableSelection` (which chooses them)
each carried their own copy of all three, byte-identical apart from whitespace. Neither module
depends on the other, so an edit to the group list in one did not reach the other - unlike
`MASK_METHODS`, whose duplication is deliberate and says so in place, this one was silent.

The column vectors are shared rather than copied per module. Every use in both modules reads them
(`GROUP_COLUMNS[group]` into a loop, a `join`, a `reduce(vcat, ...)` that builds a new vector), so
there is nothing to alias; do not start mutating them.
"""
module CovariateGroups

export JOINT_GROUP_ORDER, GROUP_COLUMNS, GROUP_FAMILY

"""The nine covariate groups, in the order the selection and the design both iterate them."""
const JOINT_GROUP_ORDER = [
    "elevation", "slope", "aspect", "t2m_c", "d2m_c",
    "u10", "v10", "sp_hpa", "ndvi",
]

"""Each group's columns in the joint panel. `aspect` is the only group wider than one column."""
const GROUP_COLUMNS = Dict(
    "elevation" => [:elevation_m],
    "slope" => [:slope_deg],
    "aspect" => [:aspect_sin, :aspect_cos],
    "t2m_c" => [:t2m_c],
    "d2m_c" => [:d2m_c],
    "u10" => [:u10],
    "v10" => [:v10],
    "sp_hpa" => [:sp_hpa],
    "ndvi" => [:ndvi],
)

"""Which data source a group came from, for the per-family reporting the selection writes."""
const GROUP_FAMILY = Dict(
    group => (group in ("elevation", "slope", "aspect") ? "dem" :
        group == "ndvi" ? "ndvi" : "era5") for group in JOINT_GROUP_ORDER
)

end # module
