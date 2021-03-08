if !(@isdefined BetadiversityHotspots)
    import Pkg; Pkg.activate(".")
    @time include("required.jl")
end

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
                    # title = "$(titlecase(outcome)) matrix Y (sorted by row and column sum)",
                    ylabel = "Sites", yticks = :none, xlabel = "Species number"
                    )

## Richness

# Get richness
richness = calculate_richness(Y, distributions[1])

# Plot richness
richness_plot = plotSDM2(richness; c=:viridis,
                        # title = "Richness ($outcome distributions)",
                        clim = (0.0, maximum(richness)),
                        colorbar_title = "Species richness"
                        )
# Plot richness quantiles
richness_qplot = plotSDM2(quantiles(richness), c=:viridis,
                         # title = "Richness quantiles ($outcome distributions)",
                         colorbar_title = "Species richness (quantiles)"
                         )

## LCBD

# Get LCBD values
lcbd = calculate_lcbd(Y, distributions[1]; transform = true, relative = true)

# Get total beta
beta_total = calculate_BDtotal(Y)

# Plot relative values
lcbdtr_plot = plotSDM2(lcbd, c=:viridis,
                      # title = "LCBD values per site ($(outcome) distributions, hellinger transformed)",
                      colorbar_title = "Relative LCBD value"
                      )
plot!([NaN], label = "BDtot = $(round(beta_total; digits = 5))",
      legend = :bottomright)

# Plot quantile scores
lcbdtr_qplot = plotSDM2(quantiles(lcbd), c=:viridis,
                       # title = "LCBD quantiles ($(outcome) distributions, hellinger transformed)",
                       colorbar_title = "Relative LCBD value (quantiles)"
                       )

## Relationship

# Plot relationship as histogram2d
rel2d_plot = histogram2d(richness, lcbd, c = :viridis, bins = 40, # title = "Relationship",
                         xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
                         xlim = (1.0, Inf), ylim = (0.0, 1.0),
                         # right_margin = 3.0mm,
                         # colorbar_titlefonthalign = :left,
                         size = (650, 400),
                         aspect_ratio = 40
                         )
vline!([median(richness)], label = :none, 
       linestyle = :dash, c = :grey)
hline!([median(lcbd)], label = :none, 
       linestyle = :dash, c = :grey)

## Export figures

# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome))"
    savefig(plot(Yobs_plot, dpi = 200),     joinpath("fig", outcome, "04-1_$(outcome)_Yobs.png"))
    savefig(plot(richness_plot, dpi = 200), joinpath("fig", outcome, "04-2_$(outcome)_richness.png"))
    savefig(plot(lcbdtr_plot, dpi = 200),   joinpath("fig", outcome, "04-3_$(outcome)_lcbd-transf.png"))
    savefig(plot(rel2d_plot, dpi = 200),    joinpath("fig", outcome, "04-4_$(outcome)_relationship2d-transf.png"))
else
    @info "Figures not saved ($(outcome))"
end

# save_quantile_figures = true # should quantile figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Quantile figures saved ($(outcome))"
    savefig(plot(richness_qplot, dpi = 200), joinpath("fig", "quantiles", "04-2_$(outcome)_richness_quantiles.png"))
    savefig(plot(lcbdtr_qplot, dpi = 200),   joinpath("fig", "quantiles", "04-3_$(outcome)_lcbd-transf_quantiles.png"))
else
    @info "Quantile figures not saved ($(outcome) richness)"
end
