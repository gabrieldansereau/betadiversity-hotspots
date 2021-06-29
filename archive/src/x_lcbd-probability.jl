@time include("src/required.jl")

outcome = "bart"

@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

@load joinpath("data", "jld2", "$(outcome)-distributions_xtras.jld2") prob_distrib lower_distrib upper_distrib

prob_distrib[1].grid
filter(!isnothing, prob_distrib[1].grid)

plot_layer(prob_distrib[1])
plot_layer(prob_distrib[1]; c=:BuPu)

# Create matrix Y
Y = Ymatrix(prob_distrib; transform=false)

# Create matrix of observed sites only
Yobs = _Yobs(Y)

# Sort Yobs by rows & columns sums
Ysort = sortslices(Yobs; dims=1, by=sum) |> y -> sortslices(y; dims=2, by=sum, rev=true)
# Heatmap of Yobs
Yobs_plot = heatmap(
    Ysort;
    # title = "$(titlecase(outcome)) matrix Y (sorted by row and column sum)",
    ylabel="Sites",
    yticks=:none,
    xlabel="Species number",
)

## Richness

# Get richness
richness = calculate_richness(Y, prob_distrib[1])
richness.grid

# Plot richness
richness_plot = plot_layer(
    richness;
    c=:viridis,
    # title = "Richness ($outcome distributions)",
    clim=(0.0, maximum(richness)),
    colorbar_title="Number of species per site",
)

## LCBD

# Get LCBD values
lcbd = calculate_lcbd(Y, prob_distrib[1]; transform=true, relative=true)

# Plot relative values
lcbdtr_plot = plot_layer(
    lcbd;
    c=:viridis,
    # title = "LCBD values per site ($(outcome) distributions, hellinger transformed)",
    colorbar_title="Relative LCBD value",
)

## Relationship

# Plot relationship as histogram2d
rel2d_plot = histogram2d(
    richness,
    lcbd;
    c=:viridis,
    bins=40, # title = "Relationship",
    xlabel="Richness",
    ylabel="LCBD",
    colorbar_title="Number of sites",
    xlim=(0.0, Inf),
    ylim=(0.0, 1.0),
)
