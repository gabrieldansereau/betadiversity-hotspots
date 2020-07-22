#### Master file
if !(@isdefined BetadiversityHotspots)
    import Pkg; Pkg.activate(".")
    @time include("required.jl")
end

## Data preparation

# save_prepdata = true # should preparation data be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

# 0a. EBD data extraction (if needed) (possibly very long)
@time run("Rscript src/00a_preparation_ebd-extraction.R")

# 0b. EBD data preparation (if needed)
@time include("00b_preparation_ebd-preparation.jl")

# 0c. Landcover data preparation (if needed)
@time include("00c_preparation_landcover.jl")

## Distributions & modelling

outcome = "raw" # desired outcome (required) # "raw", "bio", or "rf"
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
# save_data = true # should data files be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

# 1. Generate presence-absence data & single species map
@time include("01_distributions.jl")

# 2. Train Random Forest models
@time run("Rscript src/02_training_random-forests.R")

# 3. Get Random Forest predictions
@time include("03_predictions_random-forests.jl")

## Analyses

# outcome = "raw"
# save_figures = true # should figures be overwritten (optional)
# save_quantile_figures = true # should figures be overwritten (optional)

# 4. Run main analyses (Y-matrix, richness, LCBD, relationship)
@time include("04_analysis.jl")

# 5. Run subareas analyses
@time include("05_subareas.jl")

# 6. Run moving-windows analyses for LCBD values
@time include("06_moving-windows.jl")
