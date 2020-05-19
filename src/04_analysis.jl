import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# save_figures = true # should figures be overwritten (optional)
# save_quantile_figures = true # should quantile figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
    @warn "'outcome' not defined, must be either 'raw', 'bio', 'rf, or 'bart'"
elseif !(outcome in ("raw", "bio", "rf", "bart"))
    @warn "'outcome' invalid, must be either 'raw', 'bio', 'rf', or 'bart'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Load presence-absence data for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

## Y matrix (site-by-species community data table)

# Create matrix Y
Y = calculate_Y(distributions; transform = false)

# Create matrix of observed sites only
Yobs = _Yobs(Y)

# Sort Yobs by rows & columns sums
Ysort = sortslices(Yobs, dims = 1, by = sum) |> y ->
            sortslices(y, dims = 2, by = sum, rev = true)
# Heatmap of Yobs
Yobs_plot = heatmap(Ysort,
                    title = "$(titlecase(outcome)) matrix Y (sorted by row and column sum)",
                    ylabel = "Sites", yticks = :none, xlabel = "Species number", # dpi=300
                    )

## Richness

# Get richness
richness = calculate_richness(Y, distributions[1])

# Plot richness
richness_plot = plotSDM(richness; c=:viridis,
                        title = "Richness ($outcome distributions)",
                        clim=(0.0, Inf),
                        colorbar_title = "Number of species per site",
                        # dpi=300
                        )
# Plot richness quantiles
richness_qplot = plotSDM(quantiles(richness), c=:viridis,
                         title = "Richness quantiles ($outcome distributions)",
                         colorbar_title = "Number of species per site (quantiles)",
                         # dpi=300
                         )

## LCBD

# Get LCBD values
lcbd = calculate_lcbd(Y, distributions[1]; transform = true, relative = true)

# Plot relative values
lcbdtr_plot = plotSDM(lcbd, c=:viridis,
                      title = "LCBD values per site ($(outcome) distributions, hellinger transformed)",
                      colorbar_title = "LCBD value (relative to maximum)", # dpi=300
                      )

# Plot quantile scores
lcbdtr_qplot = plotSDM(quantiles(lcbd), c=:viridis,
                       title = "LCBD quantiles ($(outcome) distributions, hellinger transformed)",
                       colorbar_title = "LCBD quantile score", # dpi=300
                       )

## Relationship

# Plot relationship as histogram2d
rel2d_plot = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
                         xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
                         xlim = (1, 45), ylim = (0.0, 1.0), # dpi = 150
                         )


## Export figures

# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome)"
    savefig(Yobs_plot,     joinpath("fig", outcome, "04-1_$(outcome)_Yobs.png"))
    savefig(richness_plot, joinpath("fig", outcome, "04-2_$(outcome)_richness.png"))
    savefig(lcbdtr_plot,   joinpath("fig", outcome, "04-3_$(outcome)_lcbd-transf.png"))
    savefig(rel2d_plot,    joinpath("fig", outcome, "04-4_$(outcome)_relationship2d-transf.png"))
else
    @info "Figures not saved ($(outcome))"
end

# save_quantile_figures = true # should quantile figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Quantile figures saved ($(outcome))"
    savefig(richness_qplot, joinpath("fig", "quantiles", "04-2_$(outcome)_richness_quantiles.png"))
    savefig(lcbdtr_qplot,   joinpath("fig", "quantiles", "04-3_$(outcome)_lcbd-transf_quantiles.png"))
else
    @info "Quantile figures not saved ($(outcome) richness)"
end
