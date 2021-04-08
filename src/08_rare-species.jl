import Pkg; Pkg.activate(".")
include("required.jl")

## Conditional arguments
outcome = "raw"
# outcome = "bart"
# save_figures = true

## Occupancy
# Prep values
include("04_analysis.jl")
include("05_subareas.jl")

# Count sites with occurrences per species
getoccupancy(Ymatrix) = vec(sum(_Yobs(Ymatrix), dims = 1))
occupancy = getoccupancy(Y)
occupancy_NE = getoccupancy(Y_NE)
occupancy_SW = getoccupancy(Y_SW)
# Visualize
scatter(occupancy, formatter = :plain, legend = :none,
        xlabel = "Species ID number", ylabel = "Species occupancy")
scatter!(occupancy_NE)
scatter!(occupancy_SW)

# Occupancy percentage
getoccupancy_percentage(Ymatrix) = getoccupancy(Ymatrix) ./ size(_Yobs(Ymatrix), 1)
occupancy_percentage = getoccupancy_percentage(Y)
occupancy_percentage_NE = getoccupancy_percentage(Y_NE)
occupancy_percentage_SW = getoccupancy_percentage(Y_SW)
# Classify as rare species
# thresholds = [0.1, 0.2, 0.3, 0.4, 0.5]
# rarespecies = [occupancy_percentage .< t for t in thresholds]
getrarespecies(occupancy_percentage, threshold) = [iszero(op) ? missing : (op < threshold ? 1 : 0 ) for op in occupancy_percentage]
rarespecies = getrarespecies(occupancy_percentage, 0.4)
rarespecies_NE = getrarespecies(occupancy_percentage_NE, 0.4)
rarespecies_SW = getrarespecies(occupancy_percentage_SW, 0.4)
# Percentage of rarespecies
rarespecies_percentage = mean(rarespecies)
rarespecies_percentage_NE = mean(skipmissing(rarespecies_NE))
rarespecies_percentage_SW = mean(skipmissing(rarespecies_SW))

# LCBD-richness correlation
get_eusrr(richness, lcbd) = corspearman(collect(richness), collect(lcbd))
eusrr = get_eusrr(richness, lcbd)
eusrr_NE = get_eusrr(richness_NE, lcbd_NE)
eusrr_SW = get_eusrr(richness_SW, lcbd_SW)

# Plot EUSRR ~ rare species percentage
# Similar to Fig. 3 of Yao et al. 2021
# i_thr = findall(x -> x == 0.4, thresholds)
# scatter(rarespecies_percentage[i_thr], [eusrr],
scatter([rarespecies_percentage], [eusrr],
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = permutedims(string.(thresholds[i_thr])),
        legendtitle = "Thresholds", legendtitlefontsize = 9,
        )
hline!([0.0], style = :dash, c = :grey, label = :none)
# Even at high percentage of rare species, EUSRR can be negative

# Different thresholds
# scatter(permutedims(rarespecies_percentage), repeat([eusrr], 5),
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
scatter(permutedims(repeat([rarespecies_percentage], 3)), permutedims([eusrr, eusrr_NE, eusrr_SW]),
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
scatter(permutedims([rarespecies_percentage, rarespecies_percentage_NE, rarespecies_percentage_SW]), 
        permutedims([eusrr, eusrr_NE, eusrr_SW]),
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = ["Total" "NE" "SW"],
        legendtitle = "Regions", legendtitlefontsize = 9,
        ) |>
        x -> hline!(x, [0.0], style = :dash, c = :grey, label = :none)
# EUSRR relative whatever percentage of rare species?
savefig(joinpath("fig", outcome, "08_$(outcome)_rare-species_eusrr.png"))