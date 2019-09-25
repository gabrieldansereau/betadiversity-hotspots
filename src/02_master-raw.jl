#### Master file - Raw observations
## NOT TESTED

using Distributed
using JLD2
@time @everywhere include("src/required.jl")

# Precompile plot function
@time plot() # speeds thing up

## 0. EBD data preparation (if needed)
@time include("01_data-preparation.jl")

## 1. Generate presence-absence data & single species map
@time include("raw/01_raw_presence-absence.jl")

## 1b. Sample checklists within pixels for alternative presence-absence data
# NOT READY
# @time include("raw/01b_raw_pixel-sampling.jl")

## 2. Generate Y matrices data & heatmap visualization
@time include("raw/02_raw_Y-matrix.jl")

## 3. Run species richness analysis
@time include("raw/03_raw_richness.jl")

## 4. Run diversity analysis
@time include("raw/04_raw_community.jl")

## 5. Run LCBD analysis
@time include("raw/05_raw_lcbd.jl")

## 6. Run LCBD-richness relationship analysis
@time include("raw/06_raw_relation-lcbd-richness.jl")
