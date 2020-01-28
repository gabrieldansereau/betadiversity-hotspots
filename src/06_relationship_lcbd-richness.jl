import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load raw LCBD & richness results
outcome = "sdm"
# Load richness script
@time include("../$(outcome)/03_$(outcome)_richness.jl")
# Load LCBD script
@time include("../$(outcome)/05_$(outcome)_lcbd.jl")

# Stash results
raw = (distributions = distributions,
       Y = Y,
       Yobs = Yobs,
       Ytransf = Ytransf,
       inds_obs = inds_obs,
       inds_notobs = inds_notobs,
       richness = richness,
       LCBD = LCBD)

## Load SDM LCBD & richness results
outcome = "sdm"
# Load richness script
@time include("../$(outcome)/03_$(outcome)_richness.jl")
# Load LCBD script
@time include("../$(outcome)/05_$(outcome)_lcbd.jl")

# Stash results
sdm = (distributions = distributions,
       Y = Y,
       Yobs = Yobs,
       Ytransf = Ytransf,
       inds_obs = inds_obs,
       inds_notobs = inds_notobs,
       richness = richness,
       LCBD = LCBD)

## Richness-LCBD relationship
# Calculate relative richness (α/γ)
rel_richness = [res.richness.grid ./ (size(res.Y, 2)+1) for res in (raw,sdm)]
# Scatterplot LCBD ~ richness
relation_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[1].grid),
         markersize = 2,
         c = :skyblue,
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relation_plot, vec(rel_richness[2]), vec(sdm.LCBD[1].grid),
         markersize = 2, color=:orange, msw = 0, label = "SDM predictions")
relationtr_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD.transf.grid),
         markersize = 2,
         c =RGB(86/255, 180/255, 233/255),
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relationtr_plot, vec(rel_richness[2]), vec(sdm.LCBD[1].grid),
         markersize = 3, c = RGB(230/255,159/255,0/255), msw = 0, label = "SDM predictions")
relationdbtr_plot = scatter(vec(rel_richness[1]), vec(raw.LCBD[3].grid),
         markersize = 2,
         c =RGB(86/255, 180/255, 233/255),
         msw = 0,
         label = "Raw occurrence data",
         legend = :topright,
         xlims = (0.0, 1.0), ylims = (0.0, 1.0),
         yticks = 0.0:0.20:1.0,
         xlabel = "Species richness (\\alpha\\/\\gamma)", ylabel = "LCBD (relative to maximum)",
         grid=:none)
scatter!(relationdbtr_plot, vec(rel_richness[2]), vec(sdm.LCBD[3].grid),
         markersize = 3, c = RGB(230/255,159/255,0/255), msw = 0, label = "SDM predictions")

## Save result
#=
savefig(relation_plot, "fig/$(outcome)/06_$(outcome)_relation-lcbd-richness.png")
savefig(relationtr_plot, "fig/$(outcome)/06_$(outcome)_relation-lcbd-richness-transf.png")
=#
