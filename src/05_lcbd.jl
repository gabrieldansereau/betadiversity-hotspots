import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# outcome = "sdm" # desired outcome (required)
# save_figures = true # should figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
  @warn "'outcome' not defined, must be either 'raw' or 'sdm'"
elseif (outcome != "raw" && outcome != "sdm")
  @warn "'outcome' invalid, must be either 'raw' or 'sdm'"
else
  @info "'outcome' currently set to '$(outcome)'"
end

## Load distributions for all species
@load "data/jld2/$(outcome)-distributions.jld2" distributions
## Load matrix Y
@load "data/jld2/$(outcome)-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs

## Compute beta diversity statistics
# Load functions
include("lib/beta-div.jl")
# Compute BD statistics on distribution data
resBDobs = BD(Yobs)
# Compute BD statistics on transformed data
resBDtransf = BD(Ytransf)

# Extract LCBD values
resBD = [resBDobs, resBDtransf]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets_raw = copy(LCBDsets)
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]
LCBDsets = [LCBDsets..., LCBDsets_raw...]

## Arrange LCBD values as grid
# Create empty grids
LCBDgrids = [fill(NaN, size(distributions[1])) for LCBDi in LCBDsets]
# Fill in grids
[LCBDgrids[i][inds_obs] = LCBDsets[i] for i in 1:length(LCBDgrids)]
# Create SimpleSDMLayer with LCBD values
LCBD = SimpleSDMResponse.(LCBDgrids, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)

## Plot results
lcbd_plot1 = plotSDM(LCBD[1], c=:viridis)
heatmap!(lcbd_plot1, title = "LCBD values per site ($(outcome) distributions)",
         colorbar_title = "LCBD value (relative to maximum)", dpi=300)
lcbd_plot2 = plotSDM(LCBD[2], c=:viridis)
heatmap!(lcbd_plot2, title = "LCBD values per site ($(outcome) distributions, hellinger transformed)",
         colorbar_title = "LCBD value (relative to maximum)", dpi=300)

## Save result
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) lcbd)"
    savefig(lcbd_plot1, "fig/$(outcome)/05_$(outcome)_lcbd.pdf")
    savefig(lcbd_plot2, "fig/$(outcome)/05_$(outcome)_lcbd-transf.pdf")
else
    @info "Figures not saved ($(outcome) lcbd)"
end
