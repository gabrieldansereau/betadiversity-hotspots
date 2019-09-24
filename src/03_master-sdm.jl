#### Master file - SDM predictions
## NOT TESTED

using Distributed
using JLD2
@time @everywhere include("src/required.jl")

# Precompile plot function
@time plot() # speeds thing up

## 0. EBD data preparation (if needed)
@time include("data-preparation.jl")

## 1. Make SDM predictions for all species
@time include("sdm/sdm-predictions.jl")

## 1b. Make SDM predictions 1 species at the time across its whole range
# @time include("sdm/sdm-single-species.jl")

## 2. Generate Y matrices data & heatmap visualization
@time include("sdm/sdm-Y-matrix.jl")

## 3. Run species richness analysis
@time include("sdm/sdm-richness.jl")

## 4. Run diversity analysis
@time include("sdm/sdm-community.jl")

## 5. Run LCBD analysis
@time include("sdm/sdm-lcbd.jl")

## 6. Run LCBD-richness relationship analysis
@time include("sdm/sdm-relation-lcbd-richness.jl")
