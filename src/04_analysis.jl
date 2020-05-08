import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# save_figures = true # should figures be overwritten (optional)
# save_quantile_figures = true # should quantile figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
    @warn "'outcome' not defined, must be either 'raw' or 'sdm'"
elseif (outcome != "raw" && outcome != "sdm")
    @warn "'outcome' invalid, must be either 'raw' or 'sdm'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Load presence-absence data for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

## Y matrix (site-by-species community data table)

# Create matrix Y
Y = calculate_Y(distributions; transform = false)

# Heatmap of Y
Yobs = _Yobs(Y)
heat_y = heatmap(Yobs, title = "$(titlecase(outcome)) matrix Y (unsorted)",
                   ylabel = "Site number", xlabel = "Species number")
# Sort Y by rows & columns sums
rowsum = sum.(eachrow(Yobs))
colsum = sum.(eachcol(Yobs))
sortedrows = sortperm(rowsum)
sortedcols = sortperm(colsum, rev=true)
heat_sortrowcol = heatmap(Yobs[sortedrows, sortedcols],
                          title = "$(titlecase(outcome)) matrix Y (sorted by row and column sums)",
                          ylabel = "Site number", xlabel = "Species number", dpi=300)

## Richness

# Get richness
richness = calculate_richness(Y, distributions[1])

# Plot richness
richness_plot = plotSDM(richness, c=:viridis,
                        title = "Richness ($outcome distributions)",
                        clim=(0.0, Inf),
                        colorbar_title = "Number of species per site",
                        dpi=300)
# Plot richness quantiles
richness_qplot = plotSDM(quantiles(richness), c=:viridis,
                         title = "Richness quantiles ($outcome distributions)",
                         colorbar_title = "Number of species per site (quantiles)",
                         dpi=300)

## LCBD

# Get LCBD values
lcbd = calculate_lcbd(Y, distributions[1]; transform = true, relative = true)

# Plot relative values
lcbdtr_plot = plotSDM(lcbd, c=:viridis,
                      title = "LCBD values per site ($(outcome) distributions, hellinger transformed)",
                      colorbar_title = "LCBD value (relative to maximum)", dpi=300)

# Plot quantile scores
lcbdtr_qplot = plotSDM(quantiles(lcbd), c=:viridis,
                       title = "LCBD quantiles ($(outcome) distributions, hellinger transformed)",
                       colorbar_title = "LCBD quantile score", dpi=300)

## Relationship

# Plot relationship as histogram2d
rel2d = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
                    xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
                    xlim = (1, 45), ylim = (0.0, 1.0), dpi = 150)


## Export figures

# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome)"
    savefig(heat_sortrowcol, joinpath("fig", outcome, "02_$(outcome)_Y-rowcolsorted.png"))
    savefig(richness_plot, joinpath("fig", outcome, "03_$(outcome)_richness.png"))
    savefig(lcbdtr_plot, joinpath("fig", outcome, "05_$(outcome)_lcbd-transf.png"))
    savefig(rel2d, joinpath("fig", outcome, "06_$(outcome)_relationship2d-transf.png"))
else
    @info "Figures not saved ($(outcome))"
end

# save_quantile_figures = true # should quantile figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Quantile figures saved ($(outcome))"
    savefig(richness_qplot, joinpath("fig", "quantiles", "03_$(outcome)_richness_quantiles.png"))
    savefig(lcbdtr_qplot, joinpath("fig", "quantiles", "05_$(outcome)_lcbd-transf_quantiles.png"))
else
    @info "Quantile figures not saved ($(outcome) richness)"
end
