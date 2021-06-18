#### Main file ####
# Use this file to run the scripts all at once with conditional arguments
# Useful to update results & figures
# BUT be careful as some steps might take a while

# Load required packages & scripts
import Pkg; Pkg.activate(".")
include("required.jl")

## Part I - Data preparation ####

# Set conditional arguments used in scripts
save_prepdata = true # should preparation data be overwritten (optional)
save_figures = true # should figures be overwritten (optional)

# 0a. EBD data extraction (if needed) (possibly very long)
@time run("Rscript src/00a_preparation_ebd-extraction.R")

# 0b. EBD data preparation (if needed)
@time include("00b_preparation_ebd-preparation.jl")

# 0c. Landcover data preparation (if needed)
@time include("00c_preparation_landcover.jl")

## Part II - Distributions & modelling ####

# Set conditional arguments
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

## Part III - Main analyses ####

# Set conditional arguments
outcome = "raw" # desired outcome (required), choose either raw or bart
outcome = "bart"
save_figures = true # should figures be overwritten (optional)
save_quantile_figures = true # should figures be overwritten (optional)

# 4. Run main analyses (Y-matrix, richness, LCBD, relationship)
@time include("04_analysis.jl")

# 5. Run subareas analyses
@time include("05_subareas.jl")

# 6. Run moving-windows analyses for LCBD values
@time include("06_moving-windows.jl")

## Part IV - Additional analyses ####

# These analyses will re-run the previous scripts.
# Better to reset the conditional arguments to avoid overwriting things
outcome = nothing
save_figures = nothing
save_quantile_figures = nothing
save_data = nothing

# Set conditional arguments for additional analyses scripts
save_additional_figures = true
save_additional_data = true

# 7. Run comparison analyses (between raw data & predictions)
# (This script doesn't need to set the outcome)
@time include("07a_comparison_data.jl")
@time run(`Rscript src/07b_comparison_glm.R`) # this script might overflow memory
@time include("07c_comparison_plots.jl")
# Better to reset the outcome after this script
outcome = nothing

# 8. Run rare species analyses
# (This script needs one of the two outcomes)
outcome = "raw"
outcome = "bart"
@time include("08_rare-species.jl")
# Better to reset the outcome after this script
outcome = nothing