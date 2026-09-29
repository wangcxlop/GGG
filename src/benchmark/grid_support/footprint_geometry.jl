# ---------------------------------------------------------------------------------------------
# FY4B ground footprint
# ---------------------------------------------------------------------------------------------

"""Local east/north offset in km of `(lon, lat)` from `(lon0, lat0)`, for offsets of a few km."""
_local_km(lon, lat, lon0, lat0) = (
    (lon - lon0) * cosd(lat0) * (pi / 180) * EARTH_RADIUS_KM,
    (lat - lat0) * (pi / 180) * EARTH_RADIUS_KM,
)

"""
Ground size of one FY4B pixel at `(lon, lat)`, as
`(; along_scan_km, cross_scan_km, area_km2, nadir_ratio)`.

`RESOLUTION = 4000.0` in `FY4BPreprocessing` is the pixel size at the sub-satellite point. Away
from it the same scan-angle step subtends more ground, by a factor that grows with viewing angle.
This differentiates `latlon_to_xy` about the point to get the grid Jacobian, inverts it, and
measures the ground displacement of a one-pixel step along each grid axis. `nadir_ratio` is the
area relative to a 4 km square.

`h` is the finite-difference step in degrees: large enough that the scan angles differ in double
precision, small enough that the projection is linear across it.
"""
function pixel_footprint_km(
    lon::Real, lat::Real; sat_lon::Real=FY4BPreprocessing.DEFAULT_SAT_LON, h::Real=1e-3,
)
    blank = (; along_scan_km=NaN, cross_scan_km=NaN, area_km2=NaN, nadir_ratio=NaN)
    x_e, y_e = FY4BPreprocessing.latlon_to_xy(lat, lon + h, sat_lon)
    x_w, y_w = FY4BPreprocessing.latlon_to_xy(lat, lon - h, sat_lon)
    x_n, y_n = FY4BPreprocessing.latlon_to_xy(lat + h, lon, sat_lon)
    x_s, y_s = FY4BPreprocessing.latlon_to_xy(lat - h, lon, sat_lon)
    any(isnan, (x_e, y_e, x_w, y_w, x_n, y_n, x_s, y_s)) && return blank
    # J maps (dlon, dlat) to (dx, dy) in pixels.
    j11 = (x_e - x_w) / (2h); j12 = (x_n - x_s) / (2h)
    j21 = (y_e - y_w) / (2h); j22 = (y_n - y_s) / (2h)
    determinant = j11 * j22 - j12 * j21
    determinant == 0 && return blank
    # The inverse maps a one-pixel step back to (dlon, dlat).
    dlon_dx, dlat_dx = j22 / determinant, -j21 / determinant
    dlon_dy, dlat_dy = -j12 / determinant, j11 / determinant
    vx = _local_km(lon + dlon_dx, lat + dlat_dx, lon, lat)
    vy = _local_km(lon + dlon_dy, lat + dlat_dy, lon, lat)
    area = abs(vx[1] * vy[2] - vx[2] * vy[1])
    nadir = FY4B_NADIR_RESOLUTION_M / 1000
    return (; along_scan_km=hypot(vx...), cross_scan_km=hypot(vy...),
        area_km2=area, nadir_ratio=area / nadir^2)
end
