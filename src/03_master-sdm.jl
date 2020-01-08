#### Master file - SDM predictions
## NOT TESTED

Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

# Precompile plot function
@time plot() # speeds thing up

## 0. EBD data preparation (if needed)
@time include("01_data-preparation.jl")

## 1. Make SDM predictions for all species
@time include("sdm/01_sdm_predictions.jl")

## 1b. Make SDM predictions 1 species at the time across its whole range
# @time include("sdm/sdm_single-species.jl")

## 2. Generate Y matrices data & heatmap visualization
@time include("sdm/02_sdm_Y-matrix.jl")

## 3. Run species richness analysis
@time include("sdm/03_sdm_richness.jl")

## 4. Run diversity analysis
@time include("sdm/04_sdm_community.jl")

## 5. Run LCBD analysis
@time include("sdm/05_sdm_lcbd.jl")

## 5. Run LCBD permutations
#=
@time include("sdm/05_sdm_lcbd-permutations.jl")
=#

## 6. Run LCBD-richness relationship analysis
@time include("sdm/06_sdm_relation-lcbd-richness.jl")
