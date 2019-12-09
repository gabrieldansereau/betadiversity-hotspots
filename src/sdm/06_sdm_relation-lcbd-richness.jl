using Distributed
using JLD2
@time include("../required.jl")

## Load LCBD & richness scripts (if not already loaded)
#=
# Load richness script
@time include("03_sdm_richness.jl")
# Load LCBD script
@time include("05_sdm_lcbd.jl")
=#

## Richness-LCBD relationship
# Calculate relative richness (α/γ)
rel_richness = richness.grid ./ size(Y, 2)
# Scatterplot LCBD ~ richness
relation_plot = scatter(vec(rel_richness), vec(LCBD[1].grid),
        markersize = 1,
        color=:orange,
        label = "SDM predictions",
        legend = :bottomright,
        yticks = 0.0:0.20:1.0,
        title = "Relationship between LCBD and species richness",
        xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
        grid=:none)
relationtr_plot = scatter(vec(rel_richness), vec(LCBD[2].grid),
        markersize = 1,
        color=:orange,
        label = "SDM predictions",
        legend = :bottomright,
        yticks = 0.0:0.20:1.0,
        title = "Relationship between LCBD (hellinger transformed) and species richness",
        xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
        grid=:none)

## Save result
#=
savefig(relation_plot, "fig/sdm/06_sdm_relation-lcbd-richness.png")
savefig(relationtr_plot, "fig/sdm/06_sdm_relation-lcbd-richness-transf.png")
=#
