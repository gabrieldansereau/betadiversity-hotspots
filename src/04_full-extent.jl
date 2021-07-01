include("required.jl")

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# outcome = "bart"
# save_figures = true # should figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome) || !(outcome in ["raw", "bart"])
    @error "'outcome' must be either 'raw' or 'bart'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Load presence-absence data for all species

# Load species names
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame)
spenames = filter(:type => ==("species"), glossary).full_name

# Load distributions
distributions = [
    geotiff(
        SimpleSDMPredictor, joinpath("data", "raster", "distributions_$(outcome).tif"), i
    ) for i in eachindex(spenames)
]

## Y matrix (site-by-species community data table)

# Create matrix Y
Y = Ymatrix(distributions; transform=false)

# Create matrix of observed sites only
Yobs = Ymatrix(distributions; transform=false, observed=true)

## Richness

# Get richness
richness_full = richness(Y, distributions[1])

# Plot richness
richness_plot = plot_layer(
    richness_full;
    c=:viridis,
    colorbar_title="Species richness",
    size=(650, 400),
    dpi=200,
)

## LCBD

# Get relative LCBD values
lcbd_rel = lcbd(Y, distributions[1]; transform=true, relative=true)

# Get non-relative values
lcbd_full = lcbd(Y, distributions[1]; transform=true, relative=false)

# Get total beta
beta_full = beta_total(Y)

# Plot relative values
lcbd_rel_plot = plot_layer(
    lcbd_rel;
    c=:viridis,
    colorbar_title="Relative LCBD value",
    size=(650, 400),
    dpi=200,
)

# Plot absolute values (scaled for a nicer result)
scaling_factor = lcbd_full |> maximum |> log10 |> abs |> ceil |> Int
scaling_value = 10^scaling_factor
lcbd_resc = rescale(lcbd_full, extrema(lcbd_full) .* scaling_value)
lcbd_plot = plot_layer(
    lcbd_resc;
    c=:viridis,
    colorbar_title="LCBD value (x $(format(scaling_value, commas=true)))",
    size=(650, 400),
    dpi=200,
)

## Relationship

# Functions to add total beta value
rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])
function rectangle!(w, h, x, y, textstring, textsize)
    return plot!(
        rectangle(w, h, x, y);
        color=:white,
        legend=:none,
        annotations=(x + (w / 2), y + (h / 2), text(textstring, textsize, :center)),
    )
end

# Plot relationship as histogram2d
lmin, lmax = extrema(lcbd_resc)
lrange = lmax - lmin

rel2d_plot = histogram2d(
    richness_full,
    lcbd_resc;
    c=:viridis,
    bins=40,
    xlabel="Richness",
    ylabel="LCBD value (x 100,000)",
    colorbar_title="Number of sites",
    xlim=(1.0, Inf),
    size=(650, 400),
    dpi=200,
)
vline!([median(richness_full)]; label=:none, linestyle=:dash, c=:grey)
hline!([median(lcbd_resc)]; label=:none, linestyle=:dash, c=:grey)
rectangle!(
    16.0,
    0.15 * lrange,
    33.0,
    lmax - 0.2 * lrange,
    "BDtot = $(round(beta_full; digits=3))",
    10,
)

## Export figures

# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome))"
    savefig(richness_plot, joinpath("fig", outcome, "04_$(outcome)_richness.png"))
    savefig(lcbd_plot, joinpath("fig", outcome, "04_$(outcome)_lcbd.png"))
    savefig(rel2d_plot, joinpath("fig", outcome, "04_$(outcome)_relationship.png"))
else
    @info "Figures not saved ($(outcome))"
end
