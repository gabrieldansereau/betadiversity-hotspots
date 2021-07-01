outcome = "bart"

include(abspath("src", "04_full-extent.jl"))

## SCBD
# Prep values
inds_obs = _indsobs(Y)
Yobs = _Yobs(Y, inds_obs)
Ytransf = _Ytransf(Yobs)
BDstats = betadiv(Ytransf)

# Get SCBD values
scbd = vec(BDstats.SCBDj)

# Visualize
histogram(scbd; bins=20)
scatter(scbd)

## Species counts
# Count sites with occurrences per species
specounts = vec(sum(Yobs; dims=1))

# Visualize
histogram(specounts; bins=20)
occupancy_plot = scatter(
    specounts;
    formatter=:plain,
    legend=:none,
    xlabel="Species ID number",
    ylabel="Species occupancy",
)

## Combine stuff
scbd_plot = scatter(
    specounts,
    scbd;
    xlabel="Species occupancy",
    ylabel="SCBD value",
    formatter=:plain,
    label=:none,
)
scatter!([NaN]; label="r = $(round(cor(specounts, scbd), digits = 3))")

## Save figures
savefig(occupancy_plot, joinpath(".", "fig", "bart", "x_sp-occupancy.png"))
savefig(scbd_plot, joinpath(".", "fig", "bart", "x_sp-scbd.png"))
