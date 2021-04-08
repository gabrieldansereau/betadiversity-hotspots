import Pkg; Pkg.activate(".")
include("required.jl")

## Conditional arguments
# outcome = "raw"
outcome = "raw"
# save_figures = true

## Occupancy
# Prep values
include("04_analysis.jl")

# Count sites with occurrences per species
occupancy = vec(sum(Yobs, dims = 1))
# Visualize
occupancy_plot = scatter(occupancy, formatter = :plain, legend = :none,
                         xlabel = "Species ID number", ylabel = "Species occupancy")

# Occupancy percentage
occupancy_percentage = occupancy ./ size(Yobs, 1)
# Classify as rare species
thresholds = [0.1, 0.2, 0.3, 0.4, 0.5]
rarespecies = [occupancy_percentage .< t for t in thresholds]
# Number of rarespecies
sum.(rarespecies)
# Percentage of rarespecies
rarespecies_percentage = sum.(rarespecies)./length.(rarespecies)

# LCBD-richness correlation
eusrr = corspearman(collect(richness), collect(lcbd))

# Plot EUSRR ~ rare species percentage
# Similar to Fig. 3 of Yao et al. 2021
scatter(permutedims(rarespecies_percentage), repeat([eusrr], 5),
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.40, 1.0), ylim = (-1.0, 1.0),
        label = permutedims(string.(thresholds)),
        legendtitle = "Thresholds", legendtitlefontsize = 9,
        )
hline!([0.0], style = :dash, c = :grey, label = :none)
