# Write corrected results back to wide CSV format: columns = stations, rows = times

function write_wide(path::AbstractString, times::Vector{DateTime}, station_ids::Vector{String}, Y::Matrix{Float64})
        n_station, n_time = size(Y)
        @assert n_station == length(station_ids)
        @assert n_time == length(times)

	df = DataFrame(time = Dates.format.(times, dateformat"yyyy-mm-ddTHH:MM:SS"))
	Yt = permutedims(Y) # [ntime, nstation]
        for (j, sid) in enumerate(station_ids)
                df[!, Symbol(sid)] = Yt[:, j]
        end
        tmp_path = string(path, ".tmp-", getpid(), ".csv")
        try
                CSV.write(tmp_path, df)
                mv(tmp_path, path; force=true)
        catch e
                isfile(tmp_path) && rm(tmp_path; force=true)
                rethrow(e)
        end
end


function write_fullfit_product(
        outdir::AbstractString, product::AbstractString, times::Vector{DateTime}, ids::Vector{String},
		lonlat::Matrix{Float64}, Y_obs::Matrix{Float64}, Y_sat::Matrix{Float64}, dMat::Matrix{Float64}, best,
)
	wMat = gw_weight(dMat, Float64(best.bw); kernel=Int(best.kernel), adaptive=Bool(best.adaptive))
	R = Y_obs .- Y_sat
	Rhat = local_linear_residual_predict(
		lonlat, R, lonlat, wMat; slope_ridge=Float64(best.slope_ridge),
	)
	Y_corr = Y_sat .+ Rhat
        path = joinpath(outdir, "corr_$(product)_fullfit_insample.csv")
        try
                write_wide(path, times, ids, Y_corr)
                return path
        catch e
                stamp = Dates.format(now(), dateformat"yyyymmdd_HHMMSS")
                fallback_path = joinpath(outdir, "corr_$(product)_fullfit_insample_$(stamp).csv")
                @warn "Could not replace full-fit product file; writing timestamped fallback instead." product path fallback_path exception=(e, catch_backtrace())
                write_wide(fallback_path, times, ids, Y_corr)
                return fallback_path
        end
end


function write_validation_scope(
	outdir::AbstractString; cv_scheme::AbstractString, use_loocv_eval::Bool,
	fold_scheme::Union{Nothing,Symbol}=nothing,
	slope_ridge_candidates::Vector{Float64},
	min_scan_coverage::Float64,
	kernel_candidates::Vector{Int}=Int[],
)
	# The held-out geometry decides what the run can claim. `:random` interleaves validation
	# stations with training ones, so it supports "stations not used in fitting" but NOT
	# spatial generalisation to gauge-free ground; only `:spatial_block` holds out contiguous
	# area. Asserting the spatial claim under `:random` is a scientific overclaim, not a
	# naming slip, so the string is derived rather than fixed.
	supported_claim = if fold_scheme === :spatial_block
		"gauge-network-assisted spatial bias correction/interpolation for stations in held-out spatial blocks"
	elseif fold_scheme === :random
		"gauge-network-assisted bias correction/interpolation for stations not used in fitting, " *
			"held out at random and therefore interleaved with training stations; generalisation " *
			"to gauge-free regions is NOT supported by this split"
	else
		"gauge-network-assisted bias correction/interpolation for stations not used in fitting"
	end
	rows = [
		(key="validation_target", value="held-out station locations at matched observation/satellite timestamps"),
		(key="training_signal", value="same-timestamp observed-minus-satellite residuals from training stations"),
		(key="residual_definition", value="R = P_obs - P_sat; P_corr = P_sat + R_hat"),
		(key="local_model", value="target-centred local linear weighted ridge regression"),
		(key="local_coordinates", value="east/north offsets in kilometres, scaled by weighted RMS distance"),
		(key="slope_ridge_candidates", value=join(slope_ridge_candidates, ",")),
		(key="min_scan_coverage", value=string(min_scan_coverage)),
		(key="kernel_selection", value=length(kernel_candidates) > 1 ?
			"selected jointly with bandwidth/ridge via training-fold LOOCV, candidates=$(join(kernel_candidates, ","))" :
			"fixed by caller, not selected (kernel_candidates=$(join(kernel_candidates, ",")))"),
		(key="parameter_selection", value="coverage threshold, then LOOCV RMSE, MAE, slope_ridge, bandwidth, fixed before adaptive"),
		(key="negative_precipitation_policy", value="retain raw negative corrected values without clipping"),
		(key="validation_station_observations_used_for_fit", value="false"),
		(key="same_timestamp_training_station_observations_required", value="true"),
		(key="temporal_holdout", value="false"),
		(key="parameter_scan_use_loocv_eval", value=string(use_loocv_eval)),
		(key="cv_scheme", value=cv_scheme),
		(key="fold_scheme", value=fold_scheme === nothing ? "none" : string(fold_scheme)),
		(key="supported_claim", value=supported_claim),
		(key="unsupported_claim", value="standalone satellite correction or temporal forecast without concurrent training-station observations"),
	]
	CSV.write(joinpath(outdir, "validation_scope.csv"), DataFrame(rows))
end
