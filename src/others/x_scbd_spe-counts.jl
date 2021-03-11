outcome = "bart"

include(abspath("src", "04_analysis.jl"))

## SCBD
# Prep values
inds_obs = _indsobs(Y)
Yobs = _Yobs(Y, inds_obs)
Ytransf = _Ytransf(Yobs)
BDstats = BD(Ytransf)

# Get SCBD values
scbd = vec(BDstats.SCBDj)

# Visualize
histogram(scbd, bins = 20)
scatter(scbd)

## Species counts
# Count sites with occurrences per species
specounts = vec(sum(Yobs, dims = 1))

# Visualize
scatter(specounts)

## Combine stuff
scatter(specounts, scbd, 
        xlabel = "Species occupancy", ylabel = "SCBD value",
        formatter = :plain, label = :none)
scatter!([NaN], label = "r = $(round(cor(specounts, scbd), digits = 3))")
