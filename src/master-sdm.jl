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
@time include("sdm-predictions.jl")

## 1b. Make SDM predictions 1 species at the time across its whole range
# @time include("sdm-main-warblers.jl")

## 2. Generate Y matrices data & heatmap visualization
@time include("sdm-Y-matrix.jl")

## 3. Run species richness analysis
@time include("sdm-richness-warblers.jl")

## 4. Run diversity analysis
@time include("sdm-community-warblers.jl")

## 5. Run LCBD analysis
@time include("sdm-lcbd-warblers.jl")

## 6. Run LCBD-richness relationship analysis
@time include("sdm-relation-lcbd-richness.jl")
