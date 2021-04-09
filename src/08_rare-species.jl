import Pkg; Pkg.activate(".")
include("required.jl")

## Conditional arguments
# outcome = "raw"
outcome = "bart"
# save_figures = true

## Occupancy
# Prep values
include("04_analysis.jl")
include("05_subareas.jl")

# Count sites with occurrences per species
function getoccupancy(Ymatrix; count_sites = false)
    occupancy = vec(sum(_Yobs(Ymatrix), dims = 1))
    if !count_sites
        occupancy = occupancy ./ size(_Yobs(Ymatrix), 1)
    end
    return occupancy
end
occupancy_count = getoccupancy(Y; count_sites = true)
occupancy = getoccupancy(Y)
# occupancy_NE = getoccupancy(Y_NE)
# occupancy_SW = getoccupancy(Y_SW)

# Visualize
scatter(occupancy, formatter = :plain, legend = :none,
        xlabel = "Species ID number", ylabel = "Species occupancy")

# Classify as rare species
# Based on local occupancy
function getrarespecies(occupancy::Vector{<:AbstractFloat}, threshold::AbstractFloat)
    @assert all(x -> 0.0 <= x <= 1.0, occupancy) "occupancy must be between 0 and 1"
    rarespecies = Vector{Union{Missing, Int64}}(undef, length(occupancy))
    for (i, op) in enumerate(occupancy)
        if iszero(op)
            rarespecies[i] = missing
        else
            rarespecies[i] = op < threshold ? 1 : 0
        end
    end
    return rarespecies
end
# Based on total occupancy
function getrarespecies(occupancy::Vector{<:T}, total_occupancy::Vector{<:T}, threshold::T) where {T <: AbstractFloat}
    @assert all(x -> 0.0 <= x <= 1.0, occupancy) "occupancy must be between 0 and 1"
    @assert all(x -> 0.0 <= x <= 1.0, total_occupancy) "occupancy must be between 0 and 1"
    rarespecies = Vector{Union{Missing, Int64}}(undef, length(occupancy))
    for (i, op) in enumerate(occupancy)
        if iszero(op)
            rarespecies[i] = missing
        else
            rarespecies[i] = total_occupancy[i] < threshold ? 1 : 0
        end
    end
    return rarespecies
end
threshold = 0.4
rarespecies = getrarespecies(occupancy, threshold)
# rarespecies_NE = getrarespecies(occupancy_NE, 0.4)
# rarespecies_SW = getrarespecies(occupancy_SW, 0.4)
# rarespecies_NE_total = getrarespecies(occupancy_NE, occupancy, 0.4)
# rarespecies_SW_total = getrarespecies(occupancy_SW, occupancy, 0.4)
# Percentage of rarespecies
rarespecies_p = mean(skipmissing(rarespecies))
# rarespecies_p_NE = mean(skipmissing(rarespecies_NE))
# rarespecies_p_SW = mean(skipmissing(rarespecies_SW))
# rarespecies_p_NE_total = mean(skipmissing(rarespecies_NE_total))
# rarespecies_p_SW_total = mean(skipmissing(rarespecies_SW_total))

# Wrap as function
function get_rarespecies_p(Y, threshold)
    occupancy = getoccupancy(Y)
    rarespecies = getrarespecies(occupancy, threshold)
    rarespecies_p = mean(skipmissing(rarespecies))
    return rarespecies_p
end
get_rarespecies_p(Y, threshold)

# LCBD-richness correlation
get_eusrr(richness, lcbd) = corspearman(collect(richness), collect(lcbd))
eusrr = get_eusrr(richness, lcbd)
eusrr_NE = get_eusrr(richness_NE, lcbd_NE)
eusrr_SW = get_eusrr(richness_SW, lcbd_SW)

# Thresholds variation only
thresholds = [0.1, 0.2, 0.3, 0.4, 0.5]
[get_rarespecies_p(Y, t) for t in thresholds]
# Thresholds and region 
Ys = [Y, Y_NE, Y_SW]
rarespecies_matrix = [get_rarespecies_p(y, t) for t in thresholds, y in Ys]

# Scaling EUSRR
threshold = 0.4
rarespecies_scaling = Vector{Float64}()
eusrr_scaling = Vector{Float64}()
for sc in subarea_coords
    distributions_sc = [d[sc] for d in distributions]
    Y_sc = calculate_Y(distributions_sc)
    richness_sc = calculate_richness(Y_sc, distributions_sc[1])
    lcbd_sc = calculate_lcbd(Y_sc, distributions_sc[1])
    
    rarespecies_sc = get_rarespecies_p(Y_sc, threshold)
    eusrr_sc = get_eusrr(richness_sc, lcbd_sc)
    
    push!(rarespecies_scaling, rarespecies_sc)
    push!(eusrr_scaling, eusrr_sc)
end
rarespecies_scaling
eusrr_scaling

# Plot EUSRR ~ rare species percentage
# Similar to Fig. 3 of Yao et al. 2021
scatter([rarespecies_p], [eusrr],
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = threshold,
        legendtitle = "Thresholds", legendtitlefontsize = 9,
        )
hline!([0.0], style = :dash, c = :grey, label = :none)
# Even at high percentage of rare species, EUSRR can be negative

# Different thresholds
# scatter(permutedims(rarespecies_p), repeat([eusrr], 5),
#         xlabel = "Percentage of rare species (%)",
#         ylabel = "EUSRR",
#         xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
#         label = permutedims(string.(thresholds)),
#         legendtitle = "Thresholds", legendtitlefontsize = 9,
#         )
# hline!([0.0], style = :dash, c = :grey, label = :none)
# EUSRR doesn't vary with the threshold
# But with a lower threshold, it would be in the same area as Yao's negative EUSRR

# Subareas
scatter(permutedims(repeat([rarespecies_p], 3)), 
        [eusrr eusrr_NE eusrr_SW]),
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# EUSRR slightly more positive for subareas, but still negative
# BUT that's not quite right, how do I get the percentage of rare species in the subregion?

# Subareas
scatter(permutedims([rarespecies_p, rarespecies_p_NE, rarespecies_p_SW]), 
        permutedims([eusrr, eusrr_NE, eusrr_SW]),
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# EUSRR negative whatever subarea & percentage of rare species?
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_eusrr.png"))

# Subareas & rarity based on total occupancy
scatter(permutedims([rarespecies_p, rarespecies_p_NE_total, rarespecies_p_SW_total]), 
        permutedims([eusrr, eusrr_NE, eusrr_SW]),
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# EUSRR negative whatever subarea & percentage of rare species?
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_eusrr_total.png"))

# Subareas & thresholds
scatter(rarespecies_matrix, 
        [eusrr eusrr_NE eusrr_SW],
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.0, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# EUSRR negative whatever subarea, percentage of rare species & threshold
# Whatever the threshold, EUSRR negative
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_eusrr_thresholds.png"))

plot(thresholds, rarespecies_matrix,
     xlabel = "Rare species threshold",
     ylabel = "Rare species percentage",
     labels = ["Total" "NE" "SW"],
     legend = :bottomright)
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_thresholds.png"))

# Scaling EUSRR
scatter(rarespecies_scaling, 
        eusrr_scaling,
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.0, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# Huh ok...
show(stdout, "text/plain", rarespecies_scaling)
show(stdout, "text/plain", eusrr_scaling)
# Whatever
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_scaling.png"))


## Spatial distribution of rare species proportion
function get_site_rarespecies(Y, rarespecies, richness)
    if any(ismissing, rarespecies)
        rarespecies = replace(rarespecies, missing => 0)
    end
    Yobs = _Yobs(Y)
    Yrare = Yobs + repeat(permutedims(rarespecies), size(Yobs, 1)) .- 1.0
    Yrare_p = sum(isone, Yrare, dims = 2) ./ sum(Yobs, dims = 2)

    rarespecies_layer = copy(richness)
    inds_obs = findall(!isnothing, rarespecies_layer.grid)
    rarespecies_layer.grid[inds_obs] .= vec(Yrare_p)
    return rarespecies_layer
end
rarespecies_layer = get_site_rarespecies(Y, rarespecies, richness)
plotSDM2(rarespecies_layer, c = :viridis)
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_spatial_global.png"))

rarespecies_layer_NE = get_site_rarespecies(Y_NE, rarespecies_NE, richness_NE)
p1 = plotSDM2(rarespecies_layer_NE, c = :viridis)

rarespecies_layer_NE_total = get_site_rarespecies(Y_NE, rarespecies_NE_total, richness_NE)
p2 = plotSDM2(rarespecies_layer_NE_total, c = :viridis)

rarespecies_layer_SW = get_site_rarespecies(Y_SW, rarespecies_SW, richness_SW)
p3 = plotSDM2(rarespecies_layer_SW, c = :viridis)

rarespecies_layer_SW_total = get_site_rarespecies(Y_SW, rarespecies_SW_total, richness_SW)
p4 = plotSDM2(rarespecies_layer_SW_total, c = :viridis)

plot(p1, p2, p3, p4, dpi = 200)
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_spatial_subareas.png"))

# Rare species percentage - LCBD relationship?
histogram2d(rarespecies_layer, lcbd, c = :viridis,
            xlabel = "Percentage of rare species", ylabel = "LCBD")
histogram2d(rarespecies_layer_NE, lcbd_NE, c = :viridis,
            xlabel = "Percentage of rare species", ylabel = "LCBD")
histogram2d(rarespecies_layer_SW, lcbd_SW, c = :viridis,
            xlabel = "Percentage of rare species", ylabel = "LCBD")
# Same as scatterplots ?!
scatter(collect(rarespecies_layer), collect(lcbd))
scatter(collect(rarespecies_layer_NE), collect(lcbd_NE))
scatter(collect(rarespecies_layer_SW), collect(lcbd_SW))
scatter(collect(rarespecies_layer_NE_total), collect(lcbd_NE))
scatter(collect(rarespecies_layer_SW_total), collect(lcbd_SW))