import Pkg; Pkg.activate(".")
include("required.jl")

## Conditional arguments
# outcome = "raw"
# outcome = "bart"
# save_figures = true

# Make sure "outcome" is defined
if !(@isdefined outcome) || !(outcome in ["raw", "bart"])
    @error "'outcome' must be either 'raw' or 'bart'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

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
occupancy_NE = getoccupancy(Y_NE)
occupancy_SW = getoccupancy(Y_SW)

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
rarespecies_NE = getrarespecies(occupancy_NE, 0.4)
rarespecies_SW = getrarespecies(occupancy_SW, 0.4)
rarespecies_NE_total = getrarespecies(occupancy_NE, occupancy, 0.4)
rarespecies_SW_total = getrarespecies(occupancy_SW, occupancy, 0.4)
# Percentage of rarespecies
rarespecies_p = mean(skipmissing(rarespecies))
rarespecies_p_NE = mean(skipmissing(rarespecies_NE))
rarespecies_p_SW = mean(skipmissing(rarespecies_SW))
rarespecies_p_NE_total = mean(skipmissing(rarespecies_NE_total))
rarespecies_p_SW_total = mean(skipmissing(rarespecies_SW_total))

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

# Plot EUSRR ~ rare species percentage
# Similar to Fig. 3 of Yao et al. 2021
scatter([rarespecies_p rarespecies_p_NE rarespecies_p_SW], 
        [eusrr eusrr_NE eusrr_SW],
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# EUSRR negative whatever subarea & percentage of rare species?

# Subareas & rarity based on total occupancy
scatter([rarespecies_p rarespecies_p_NE_total rarespecies_p_SW_total], 
        [eusrr eusrr_NE eusrr_SW],
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# EUSRR negative whatever subarea & percentage of rare species?

## Effect of thresholds
# Thresholds variation only
thresholds = [0.1, 0.2, 0.3, 0.4, 0.5]
[get_rarespecies_p(Y, t) for t in thresholds]
# Thresholds and region 
Ys = [Y, Y_NE, Y_SW]
rarespecies_matrix = [get_rarespecies_p(y, t) for t in thresholds, y in Ys]

# Subareas & thresholds
p_eusrr = scatter(rarespecies_matrix, 
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
if (@isdefined save_figures) && save_figures == true
    savefig(p_eusrr, joinpath("fig", outcome, "08_$(outcome)_rare-species_eusrr_thresholds.png"))
end

p_thresholds = plot(thresholds, rarespecies_matrix,
     xlabel = "Rare species threshold",
     ylabel = "Rare species percentage",
     labels = ["Total" "NE" "SW"],
     legend = :bottomright)

if (@isdefined save_figures) && save_figures == true
    savefig(p_thresholds, joinpath("fig", outcome, "08_$(outcome)_rare-species_thresholds.png"))
end

## Effect of scaling
# Scaling EUSRR
rarespecies_scaling = Vector{Float64}()
eusrr_scaling = Vector{Float64}()
for sc in subarea_coords
    distributions_sc = [d[sc] for d in distributions]
    Y_sc = calculate_Y(distributions_sc)
    richness_sc = calculate_richness(Y_sc, distributions_sc[1])
    lcbd_sc = calculate_lcbd(Y_sc, distributions_sc[1]; relative = false)
    
    rarespecies_sc = get_rarespecies_p(Y_sc, threshold)
    eusrr_sc = get_eusrr(richness_sc, lcbd_sc)
    
    push!(rarespecies_scaling, rarespecies_sc)
    push!(eusrr_scaling, eusrr_sc)
end
rarespecies_scaling
eusrr_scaling

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

rarespecies_layer_NE = get_site_rarespecies(Y_NE, rarespecies_NE, richness_NE)
p1 = plotSDM2(rarespecies_layer_NE, c = :viridis)

rarespecies_layer_NE_total = get_site_rarespecies(Y_NE, rarespecies_NE_total, richness_NE)
p2 = plotSDM2(rarespecies_layer_NE_total, c = :viridis)

rarespecies_layer_SW = get_site_rarespecies(Y_SW, rarespecies_SW, richness_SW)
p3 = plotSDM2(rarespecies_layer_SW, c = :viridis)

rarespecies_layer_SW_total = get_site_rarespecies(Y_SW, rarespecies_SW_total, richness_SW)
p4 = plotSDM2(rarespecies_layer_SW_total, c = :viridis)

p_subareas = plot(p1, p2, p3, p4, dpi = 200)
if (@isdefined save_figures) && save_figures == true
    savefig(p_subareas, joinpath("fig", outcome, "08_$(outcome)_rare-species_spatial_subareas.png"))
end

## Ascending & descending parts
using StatsPlots

# Check relationship plots
combined_plot

function ascending_plots(richness, lcbd, rarespecies)
    # Choose threshold
    ascending_threshold = minimum(lcbd)
    min_indx = findall(x -> x == ascending_threshold, lcbd.grid)
    abs_min = median(richness.grid[min_indx])
    
    binlayer = replace(richness, Pair.(unique(richness), unique(richness) .> abs_min)...)
    lplot = plotSDM2(binlayer, c = cgrad(:PuOr, rev = true),
                     colorbar = :none)
    
    bin_names = map(x -> isone(x) ? "Ascending" : "Descending", collect(binlayer))
    isascending = isequal.(bin_names, "Ascending")

    rare_values = collect(rarespecies)

    dplot = density(rare_values[isascending], c = :PuOr, label = "Ascending")
    density!(rare_values[.!isascending], c = cgrad(:PuOr, rev = true), label = "Descending")
    plot!(xlabel = "Rare species percentage", ylabel = "Density", bottommargin = 4.0mm)
    
    l = @layout [a{0.5w} b{0.45w}]
    plot(lplot, dplot, size = (900, 300), layout = l)
end
p_asc1 = ascending_plots(richness_NE, lcbd_NE, rarespecies_layer_NE)
p_asc2 = ascending_plots(richness_SW, lcbd_SW, rarespecies_layer_SW)
asc_plots = plot(
    p_asc1, p_asc2, 
    layout = (2,1), size = (900, 660),
    title = ["a) Northeast subregion" "" "b) Southwest subregion" ""],
    titleloc = :left,
)
if (@isdefined save_figures) && save_figures == true
    savefig(asc_plots, joinpath("fig", outcome, "08_$(outcome)_rare-species_ascending_plots.png"))
end

