import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Conditional arguments
outcome = "raw" # desired outcome

## Load presence-absence data for all species
@load "data/jld2/$(outcome)-distributions.jld2" distributions spenames speindex
