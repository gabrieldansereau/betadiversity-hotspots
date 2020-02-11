#### Master file - Raw observations
import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

# Precompile plot function
@time plot() # speeds thing up

#### Data preparation

## 0a. EBD data preparation (if needed)
# save_prepdata = true # should prepared data be overwritten (optional)
@time include("00a_data_ebd-preparation.jl")

## 0b. Landcover data preparation (if needed)
# save_figures = true # should figures be overwritten (optional)
@time include("00b_data_landcover-copernicus.jl")

#### Analyses

outcome = "raw" # desired outcome (required)
# outcome = "sdm" # desired outcome (required)
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
# save_data = true # should data files be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)
# save_relfigures = true # should relationship figures be overwritten (optional)

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

## 7. Run endemism analysis
@time incluce("07_endemism.jl")
