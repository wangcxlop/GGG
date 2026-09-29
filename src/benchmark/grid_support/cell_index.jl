# ---------------------------------------------------------------------------------------------
# Which cell a gauge reads
# ---------------------------------------------------------------------------------------------

"""
The 1-based FY4B disk pixel `FY4BPreprocessing.extract_precipitation` reads for a point, as
`(ix, iy, x, y)` where `(x, y)` is the fractional 0-based position it rounded.

`(0, 0)` marks the `NaN` that `latlon_to_xy` returns, which is rarer than it looks: the limb test
in `latlon_to_scan_angles` measures the angle *at the satellite*, and no point on Earth subtends
more than about 8.7 degrees from geostationary orbit, so it never fires and the far hemisphere
folds onto valid pixel indices instead of being rejected. What actually guards the disk is the
bounds check in `extract_precipitation`, repeated here wherever this module samples. Every gauge
in this study sits well inside the disk, so nothing below depends on the difference.
"""
function fy4b_pixel_index(lon::Real, lat::Real; sat_lon::Real=FY4BPreprocessing.DEFAULT_SAT_LON)
    x, y = FY4BPreprocessing.latlon_to_xy(lat, lon, sat_lon)
    (isnan(x) || isnan(y)) && return (0, 0, x, y)
    return (round(Int, x) + 1, round(Int, y) + 1, x, y)
end

"""
The index of the `step_deg` lat/lon cell containing a point, as `(ilon, ilat)`.

The grid is anchored at (-180, -90), which is how IMERG and GSMaP are gridded and therefore what
Earth Engine samples at `scale: 11132`. That anchoring is an assumption about the export, not
something this repository can read off the data - `identical_series_table` is what tests it.
"""
latlon_cell_index(lon::Real, lat::Real; step_deg::Real=COARSE_CELL_DEG) =
    (floor(Int, (Float64(lon) + 180) / step_deg), floor(Int, (Float64(lat) + 90) / step_deg))

"""Centre of the `step_deg` cell `(ilon, ilat)`, as `(lon, lat)`."""
latlon_cell_center(ilon::Integer, ilat::Integer; step_deg::Real=COARSE_CELL_DEG) =
    ((ilon + 0.5) * step_deg - 180, (ilat + 0.5) * step_deg - 90)
