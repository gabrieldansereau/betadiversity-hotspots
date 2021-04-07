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
rarespecies = occupancy_percentage .< 0.4
# Number of rarespecies
sum(rarespecies)
# Percentage of rarespecies
rarespecies_percentage = sum(rarespecies)/length(rarespecies)

# LCBD-richness correlation
eusrr = corspearman(collect(richness), collect(lcbd))

# Plot EUSRR ~ rare species percentage
# Similar to Fig. 3 of Yao et al. 2021
scatter([rarespecies_percentage], [eusrr],
        xlabel = "Percentage of rare species (%)",
        ylabel = "EUSRR",
        xlim = (0.50, 1.0), ylim = (-1.0, 1.0),
        label = :none)
hline!([0.0], style = :dash, c = :grey, label = :none)
