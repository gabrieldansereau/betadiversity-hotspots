import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

outcome = "raw"

## Load LCBD & richness scripts (if not already loaded)
#=
# Load richness script
@time include("03_$(outcome)_richness.jl")
# Load LCBD script
@time include("05_$(outcome)_lcbd.jl")
=#

## Richness-LCBD relationship
# Calculate relative richness (α/γ)
rel_richness = richness.grid ./ (size(Y, 2)+1)
# Scatterplot LCBD ~ richness
relation_plot = scatter(vec(rel_richness), vec(LCBD[1].grid),
        markersize = 1,
        color=:skyblue,
        label = "$(outcome) occurrence data",
        legend = :bottomright,
        yticks = 0.0:0.20:1.0,
        title = "Relationship between LCBD and species richness",
        xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
        grid=:none)
relationtr_plot = scatter(vec(rel_richness), vec(LCBD[2].grid),
        markersize = 1,
        color=:skyblue,
        label = "$(outcome) occurrence data",
        legend = :bottomright,
        yticks = 0.0:0.20:1.0,
        title = "Relationship between LCBD (hellinger transformed) and species richness",
        xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
        grid=:none)

## Save result
#=
savefig(relation_plot, "fig/$(outcome)/06_$(outcome)_relation-lcbd-richness.png")
savefig(relationtr_plot, "fig/$(outcome)/06_$(outcome)_relation-lcbd-richness-transf.png")
=#
