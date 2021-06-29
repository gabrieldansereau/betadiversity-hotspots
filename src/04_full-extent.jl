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
Yobs = _Yobs(Y)

## Richness

# Get richness
richness = calculate_richness(Y, distributions[1])

# Plot richness
richness_plot = plot_layer(
    richness;
    c=:viridis,
    colorbar_title="Species richness",
    size=(650, 400),
    dpi=200,
)

## LCBD

# Get LCBD values
lcbd_rel = calculate_lcbd(Y, distributions[1]; transform=true, relative=true)

# Get non-relative values
lcbd_abs = calculate_lcbd(Y, distributions[1]; transform=true, relative=false)
round.(Float64.(extrema(lcbd_abs)); sigdigits=4)
lcbd = lcbd_abs

# Get total beta
beta_total = calculate_BDtotal(Y)

# Plot relative values
lcbdtr_plot = plot_layer(
    lcbd_rel;
    c=:viridis,
    colorbar_title="Relative LCBD value",
    size=(650, 400),
    dpi=200,
)
# Plot absolute values
scaling_factor = lcbd |> maximum |> log10 |> abs |> ceil |> Int
scaling_value = 10^scaling_factor
lcbd_resc = rescale(lcbd, extrema(lcbd) .* scaling_value)
lcbdtr_plot = plot_layer(
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
    richness,
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
vline!([median(richness)]; label=:none, linestyle=:dash, c=:grey)
hline!([median(lcbd_resc)]; label=:none, linestyle=:dash, c=:grey)
rectangle!(
    16.0,
    0.15 * lrange,
    33.0,
    lmax - 0.2 * lrange,
    "BDtot = $(round(beta_total; digits=3))",
    10,
)

## Export figures

# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome))"
    savefig(richness_plot, joinpath("fig", outcome, "04_$(outcome)_richness.png"))
    savefig(lcbdtr_plot, joinpath("fig", outcome, "04_$(outcome)_lcbd.png"))
    savefig(rel2d_plot, joinpath("fig", outcome, "04_$(outcome)_relationship.png"))
else
    @info "Figures not saved ($(outcome))"
end
