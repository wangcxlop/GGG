using CSV

function validate_processed(root::AbstractString)
    paths = filter(path -> endswith(lowercase(path), ".csv"), readdir(root; join=true))
    types = Dict(
        :relative_humidity => Float64,
        :wind_speed => Float64,
        :sp_hpa => Float64,
        :sp_anomaly_hpa => Float64,
        :month_bjt => String,
    )

    min_rh, max_rh = Inf, -Inf
    min_wind = Inf
    min_sp, max_sp = Inf, -Inf
    max_month_abs_mean = 0.0
    bad_numeric_rows = 0
    total_rows = 0

    for path in paths
        anomaly_sums = Dict{String,Float64}()
        anomaly_counts = Dict{String,Int}()

        for row in CSV.Rows(path; types=types, reusebuffer=true)
            rh = row.relative_humidity
            wind = row.wind_speed
            sp = row.sp_hpa
            anomaly = row.sp_anomaly_hpa
            month = row.month_bjt
            total_rows += 1

            valid = all(isfinite, (rh, wind, sp, anomaly)) && 0 <= rh <= 100 && wind >= 0 && sp > 0
            bad_numeric_rows += !valid
            min_rh, max_rh = min(min_rh, rh), max(max_rh, rh)
            min_wind = min(min_wind, wind)
            min_sp, max_sp = min(min_sp, sp), max(max_sp, sp)
            anomaly_sums[month] = get(anomaly_sums, month, 0.0) + anomaly
            anomaly_counts[month] = get(anomaly_counts, month, 0) + 1
        end

        for month in keys(anomaly_sums)
            month_mean = anomaly_sums[month] / anomaly_counts[month]
            max_month_abs_mean = max(max_month_abs_mean, abs(month_mean))
        end
    end

    println("files: ", length(paths))
    println("rows: ", total_rows)
    println("bad numeric rows: ", bad_numeric_rows)
    println("relative humidity range: ", min_rh, " to ", max_rh)
    println("minimum wind speed: ", min_wind)
    println("surface pressure range (hPa): ", min_sp, " to ", max_sp)
    println("maximum absolute station-month pressure anomaly mean (hPa): ", max_month_abs_mean)

    length(paths) == 237 || error("Expected 237 station files, found $(length(paths))")
    total_rows == 237 * 26_304 || error("Unexpected total row count: $total_rows")
    bad_numeric_rows == 0 || error("Found $bad_numeric_rows invalid numeric rows")
end

root = length(ARGS) == 1 ? ARGS[1] : joinpath("data", "processed", "era5_land", "stations_2022_2024")
validate_processed(root)
