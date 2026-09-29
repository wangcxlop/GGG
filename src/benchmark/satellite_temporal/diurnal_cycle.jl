"""
First-harmonic fit `mean + amplitude * cos(2π(h - phase_hour)/24)` of a diurnal curve by least squares,
so hours with no data can simply be left out.
"""
function harmonic_phase(hours::AbstractVector{<:Real}, values::AbstractVector{<:Real})
    length(hours) == length(values) || throw(DimensionMismatch("hours and values must have the same length"))
    length(hours) >= 3 || throw(ArgumentError("need at least three hours"))
    omega = 2π .* Float64.(hours) ./ 24
    beta = [ones(length(omega)) cos.(omega) sin.(omega)] \ Float64.(values)
    return (; mean=beta[1], amplitude=hypot(beta[2], beta[3]), phase_hour=mod(atan(beta[3], beta[2]) * 24 / 2π, 24))
end

"""`a - b` on the 24-hour clock, in `[-12, 12)`."""
circular_hour_difference(a::Real, b::Real) = mod(a - b + 12, 24) - 12

"""
Composite diurnal cycles per sample x source x season x region x hour of day (hour-ending BJT).

Over the sample's cells: `mean_amount` (every hour, dry ones included), `wet_freq` (share of hours
`>= 0.1`), `wet_intensity` (mean over wet hours) and `freq_<class>` for each rain class. `sources` is
a vector of `name => matrix`, the gauge first by convention.
"""
function diurnal_cycle_table(ctx, sources::AbstractVector{<:Pair}, mask::AbstractMatrix; sample::AbstractString)
    G = length(ctx.regions) + 1
    rows = NamedTuple[]
    for (source, Y) in sources
        _check_inputs(ctx, Y, mask)
        n, n_wet = zeros(Int, 5, G, 24), zeros(Int, 5, G, 24)
        amount, wet_amount = zeros(5, G, 24), zeros(5, G, 24)
        class_n = zeros(Int, 5, G, 24, 5)
        for j in axes(Y, 2)
            s, h = ctx.column_season[j] + 1, ctx.column_hour[j] + 1
            for i in axes(Y, 1)
                mask[i, j] || continue
                v = Y[i, j]
                g = ctx.station_region_index[i] + 1
                c = intensity_class(v)
                for (a, b) in ((1, 1), (1, g), (s, 1), (s, g))
                    n[a, b, h] += 1
                    amount[a, b, h] += v
                    if c >= 1
                        n_wet[a, b, h] += 1
                        wet_amount[a, b, h] += v
                        class_n[a, b, h, c] += 1
                    end
                end
            end
        end
        for a in 1:5, b in 1:G, h in 1:24
            count_ = n[a, b, h]
            pairs = Pair{Symbol,Any}[
                :sample => sample, :source => String(source), :season => SEASON_LABELS[a],
                :region => _region_labels(ctx)[b], :hour => h - 1, :n => count_,
                :mean_amount => count_ > 0 ? amount[a, b, h] / count_ : NaN,
                :wet_freq => count_ > 0 ? n_wet[a, b, h] / count_ : NaN,
                :wet_intensity => n_wet[a, b, h] > 0 ? wet_amount[a, b, h] / n_wet[a, b, h] : NaN,
            ]
            for c in 1:5
                push!(pairs, Symbol("freq_", INTENSITY_CLASSES[c + 1]) => count_ > 0 ? class_n[a, b, h, c] / count_ : NaN)
            end
            push!(rows, _row(pairs))
        end
    end
    return DataFrame(rows)
end

"""
Each product's diurnal cycle against the gauge's, per sample x season x region.

Over the hours with data: r of the amount and frequency curves, peak hours and their circular
difference, first-harmonic phases and their difference, the relative amplitude (amplitude / mean) of
each and its ratio, and the ratio of mean amounts. `NaN` below 12 hours or 100 cells an hour.
"""
function diurnal_summary(diurnal::DataFrame, products::AbstractVector{<:AbstractString}; gauge::AbstractString="Gauge")
    rows = NamedTuple[]
    for group in groupby(diurnal, [:sample, :season, :region]; sort=true)
        reference = sort(group[group.source .== gauge, :], :hour)
        for product in products
            estimate = sort(group[group.source .== product, :], :hour)
            nrow(estimate) == nrow(reference) || continue
            keep = findall((reference.n .>= 100) .& (estimate.n .>= 100))
            hours = reference.hour[keep]
            enough = length(keep) >= 12 && sum(reference.mean_amount[keep]) > 0 && sum(estimate.mean_amount[keep]) > 0
            fit(values) = enough ? harmonic_phase(hours, values) : (; mean=NaN, amplitude=NaN, phase_hour=NaN)
            obs_amount, sat_amount = fit(reference.mean_amount[keep]), fit(estimate.mean_amount[keep])
            obs_freq, sat_freq = fit(reference.wet_freq[keep]), fit(estimate.wet_freq[keep])
            peak(values) = enough ? Float64(hours[argmax(values)]) : NaN
            obs_peak, sat_peak = peak(reference.mean_amount[keep]), peak(estimate.mean_amount[keep])
            push!(rows, (;
                sample=group.sample[1], product, season=group.season[1], region=group.region[1], n_hours=length(keep),
                r_amount=enough ? _finite_cor(reference.mean_amount[keep], estimate.mean_amount[keep]) : NaN,
                r_freq=enough ? _finite_cor(reference.wet_freq[keep], estimate.wet_freq[keep]) : NaN,
                peak_hour_gauge=obs_peak, peak_hour_sat=sat_peak,
                peak_hour_diff_h=enough ? circular_hour_difference(sat_peak, obs_peak) : NaN,
                phase_hour_gauge=obs_amount.phase_hour, phase_hour_sat=sat_amount.phase_hour,
                phase_diff_h=enough ? circular_hour_difference(sat_amount.phase_hour, obs_amount.phase_hour) : NaN,
                freq_phase_diff_h=enough ? circular_hour_difference(sat_freq.phase_hour, obs_freq.phase_hour) : NaN,
                rel_amplitude_gauge=obs_amount.amplitude / obs_amount.mean,
                rel_amplitude_sat=sat_amount.amplitude / sat_amount.mean,
                rel_amplitude_ratio=(sat_amount.amplitude / sat_amount.mean) / (obs_amount.amplitude / obs_amount.mean),
                mean_amount_ratio=obs_amount.mean > 0 ? sat_amount.mean / obs_amount.mean : NaN,
            ))
        end
    end
    return DataFrame(rows)
end
