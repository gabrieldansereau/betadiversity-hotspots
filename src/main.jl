#### Main file ####
# Use this file to run the scripts all at once with conditional arguments
# Useful to update results & figures
# BUT be careful as some steps might take a while

# Load required packages & scripts
import Pkg; Pkg.activate(".")
include("required.jl")

## Data preparation

save_prepdata = true # should preparation data be overwritten (optional)
save_figures = true # should figures be overwritten (optional)

# 0a. EBD data extraction (if needed) (possibly very long)
@time run("Rscript src/00a_preparation_ebd-extraction.R")

# 0b. EBD data preparation (if needed)
@time include("00b_preparation_ebd-preparation.jl")

# 0c. Landcover data preparation (if needed)
@time include("00c_preparation_landcover.jl")

## Distributions & modelling

outcome = "raw" # desired outcome (required) # "raw", "bio", or "rf"
create_distributions = true # should distributions be computed (optional, loaded otherwise)
save_data = true # should data files be overwritten (optional)
save_figures = true # should figures be overwritten (optional)

# 1. Generate presence-absence data & single species map
@time include("01_distributions.jl")

# 2. Train BART models
@time run("Rscript src/02_training_bart.R")

# 3. Get Random Forest predictions
@time include("03_predictions_bart.jl")

## Analyses

outcome = "raw"
outcome = "bart"
save_figures = true # should figures be overwritten (optional)
save_quantile_figures = true # should figures be overwritten (optional)

# 4. Run main analyses (Y-matrix, richness, LCBD, relationship)
@time include("04_analysis.jl")

# 5. Run subareas analyses
@time include("05_subareas.jl")

# 6. Run moving-windows analyses for LCBD values
@time include("06_moving-windows.jl")

# 7. Run comparison analyses (between raw data & predictions)
@time include("07a_comparison_data.jl")
@time run(`Rscript src/07b_comparison_glm.R`)
@time include("07c_comparison_plots.jl")

# 8. Run rare species analyses
@time include("08_rare-species.jl")