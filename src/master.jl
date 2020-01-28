#### Master file - Raw observations
## NOT TESTED

import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

# Precompile plot function
@time plot() # speeds thing up

## 0a. EBD data preparation (if needed)
@time include("00a_data_ebd-preparation.jl")

## 0b. Landcover data preparation (if needed)
@time include("00a_data_ebd-preparation.jl")

####

outcome = "raw"

## 1. Generate presence-absence data & single species map
@time include("01_distributions.jl")

## 2. Generate Y matrices data & heatmap visualization
@time include("02_Y-matrix.jl")

## 3. Run species richness analysis
@time include("03_richness.jl")

## 4. Run evenness analysis
@time include("04_evenness.jl")

## 5. Run LCBD analysis
@time include("05_lcbd.jl")

## 6. Run LCBD-richness relationship analysis
@time include("06_relationship_lcbd-richness.jl")
