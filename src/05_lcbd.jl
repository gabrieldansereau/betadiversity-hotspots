import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

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
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions
## Load matrix Y
@load joinpath("data", "jld2", "$(outcome)-Y-matrices.jld2") Y

## Compute beta diversity statistics
lcbd = calculate_lcbd(Y, distributions[1])

## Plot results
# Relative values
lcbdtr_plot = plotSDM(lcbd, c=:viridis,
                      title = "LCBD values per site ($(outcome) distributions, hellinger transformed)",
                      colorbar_title = "LCBD value (relative to maximum)", dpi=300)
# Quantile scores
lcbdtr_qplot = plotSDM(quantiles(lcbd), c=:viridis,
                       title = "LCBD quantiles ($(outcome) distributions, hellinger transformed)",
                       colorbar_title = "LCBD quantile score", dpi=300)

## Save result
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) lcbd)"
    savefig(lcbdtr_plot, joinpath("fig", outcome, "05_$(outcome)_lcbd-transf.png"))
else
    @info "Figures not saved ($(outcome) lcbd)"
end
# Quantile figures
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) lcbd)"
    savefig(lcbdtr_qplot, joinpath("fig", "quantiles", "05_$(outcome)_lcbd-transf_quantiles.png"))
else
    @info "Figures not saved ($(outcome) lcbd)"
end
