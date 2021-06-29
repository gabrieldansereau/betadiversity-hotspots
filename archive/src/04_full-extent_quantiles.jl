include("../../src/required.jl")
include("lib/quantiles.jl")

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# outcome = "bart"
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
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame)
spenames = filter(:type => ==("species"), glossary).full_name
distributions = [
    geotiff(
        SimpleSDMPredictor, joinpath("data", "proc", "distributions_$(outcome).tif"), i
    ) for i in eachindex(spenames)
]

## Y matrix (site-by-species community data table)

# Create matrix Y
Y = calculate_Y(distributions; transform=false)

# Create matrix of observed sites only
Yobs = _Yobs(Y)

## Richness

# Get richness
richness = calculate_richness(Y, distributions[1])

# Plot richness quantiles
richness_qplot = plot_layer(
    quantiles(richness);
    c=:viridis,
    # title = "Richness quantiles ($outcome distributions)",
    colorbar_title="Species richness (quantiles)",
    size=(650, 400),
)

## LCBD

# Get LCBD values
lcbd_rel = calculate_lcbd(Y, distributions[1]; transform=true, relative=true)

# Plot quantile scores
lcbdtr_qplot = plot_layer(
    quantiles(lcbd_rel);
    c=:viridis,
    # title = "LCBD quantiles ($(outcome) distributions, hellinger transformed)",
    colorbar_title="LCBD value (quantiles)",
    size=(650, 400),
)

## Export figures

# save_quantile_figures = true # should quantile figures be overwritten (optional)
if (@isdefined save_quantile_figures) && save_quantile_figures == true
    @info "Quantile figures saved ($(outcome))"
    savefig(
        plot(richness_qplot; dpi=200),
        joinpath("archive", "fig", "quantiles", "04_$(outcome)_richness_quantiles.png"),
    )
    savefig(
        plot(lcbdtr_qplot; dpi=200),
        joinpath("archive", "fig", "quantiles", "04_$(outcome)_lcbd.png"),
    )
else
    @info "Quantile figures not saved ($(outcome) richness)"
end
