#### Master file - Raw observations
## NOT TESTED

using Distributed
using JLD2
@time @everywhere include("src/required.jl")

# Precompile plot function
@time plot() # speeds thing up

## 0. EBD data preparation (if needed)
@time include("data-preparation.jl")

## 1. Generate presence-absence data & single species map
@time include("raw/raw-presence-absence.jl")

## 1b. Sample checklists within pixels for alternative presence-absence data
# NOT READY
# @time include("raw/raw-pixel-sampling.jl")

## 2. Generate Y matrices data & heatmap visualization
@time include("raw/raw-Y-matrix.jl")

## 3. Run species richness analysis
@time include("raw/raw-richness.jl")

## 4. Run diversity analysis
@time include("raw/raw-community.jl")

## 5. Run LCBD analysis
@time include("raw/raw-lcbd.jl")

## 6. Run LCBD-richness relationship analysis
@time include("raw/raw-relation-lcbd-richness.jl")
