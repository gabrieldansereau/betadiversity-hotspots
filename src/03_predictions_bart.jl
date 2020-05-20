import Pkg
Pkg.activate(".")
using RCall
begin
    R"""
    library(embarcadero)
    library(tidyverse)
    library(viridis)
    library(furrr)
    plan(multiprocess)
    """
end
using Distributed
@time include("required.jl")

## Conditional arguments
# save_data = true

## Perform BARTs
begin
    R"""
    ## 1. Load data ####

    # Load data
    spa_full <- read_tsv("data/proc/distributions_spa_full.csv")
    env_full <- read_tsv("data/proc/distributions_env_full.csv")
    spe      <- read_tsv("data/proc/distributions_spe_full.csv") 
    
    # Load QC data (optional)
    spa_qc <- read_tsv("data/proc/distributions_spa_qc.csv")
    
    # Prepare data
    subset_qc <- TRUE # subset to QC data (optional)
    source("src/02_training_data-preparation.R")

    # Select fewer variables
    xnames <- c("lat", "lon", "wc1", "wc12", paste0("lc", c(1:5, 7:10))) # lc6 always zero
    vars_df <- vars_full[,xnames]
 
    ## 2. Create layers ####

    # Convert variables to layers
    df_to_layer <- function(x, lons, lats){
        mat <- matrix(data = x, nrow = uniqueN(lats), ncol = uniqueN(lons))
        layer <- raster(
            mat[nrow(mat):1,],
            xmn=min(lons), xmx=max(lons), 
            ymn=min(lats), ymx=max(lats)
         )
        return(layer)
    }
    vars_layers <- map(
        vars_df, 
        function(x) df_to_layer(x, lons = vars_df$lon, lats = vars_df$lat)
    )
    # Stack raster layers
    (vars_stack <- stack(vars_layers, names = xnames))

    ## 3. Models ####
    
    # Load models
    load("data/proc/bart_models_qc.RData")

    # Quantile Predictions
    predictions <- future_map(
        sdms,
        function(x) predict(
            object = x, 
            x.layers = vars_stack,
            quantiles = c(0.025, 0.975),
            splitby = 20
        )
    ) # ~ 8 min., ~ 1 min in parallel

    # Predictions
    pred_df <- predictions %>% 
        map(~ .$layer.1) %>% 
        stack() %>% 
        as.data.frame(xy = TRUE) %>% 
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    pred_df
    # Lower quantiles
    lower_df <- predictions %>% 
        map(~ .$layer.2) %>% 
        stack() %>% 
        as.data.frame(xy = TRUE) %>% 
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    lower_df
    # Upper quantiles
    upper_df <- predictions %>% 
        map(~ .$layer.3) %>%
        stack() %>% 
        as.data.frame(xy = TRUE) %>%  
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    upper_df

    # Extract summary statistics
    # Inner calls
    source("src/lib/R/bart.R")
    summaries <-  map(sdms, summary_inner)

    # Organize as tibble
    results <- tibble(
        spe = names(sdms),
        auc = map_dbl(summaries, function(x) x$auc),
        threshold = map_dbl(summaries, function(x) x$threshold),
        tss = map_dbl(summaries, function(x) x$tss),
        type_I = map_dbl(summaries, function(x) x$type_I),
        type_II = map_dbl(summaries, function(x) x$type_II)
    )
    results

    # Presence-absence dataframe
    pres_df <- map2_df(
        pred_df, results$threshold, 
        function(pred, thresh) ifelse(pred > thresh, 1, 0) 
    )
    pres_df

    """
end
@rget pred_df lower_df upper_df pres_df results

## Create Y matrices

# Get matrix Y
Y = replace(Array(pres_df), missing => NaN)
Yprob  = replace(Array(pred_df),  missing => NaN)
Ylower = replace(Array(lower_df), missing => NaN)
Yupper = replace(Array(upper_df), missing => NaN)
# Set values to NaN if no species present
inds_zeros = findall(map(x -> all(iszero.(x)), eachrow(Y)))
Y[inds_zeros,:] .= NaN

## Create distributions

# Load raw distributions (for grid size)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
raw_distributions = distributions
# Cut to Quebec coordinates (optional)
coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
raw_distributions = [d[coords_qc] for d in raw_distributions]
# Get layer dimensions & limits
dims = size(raw_distributions[1].grid)
lims = (left = raw_distributions[1].left, right = raw_distributions[1].right,
        bottom = raw_distributions[1].bottom, top = raw_distributions[1].top)

# Create distribution layers
layers = []
for Y in (Y, Yprob, Ylower, Yupper)
    Ydistrib = replace(Y, 0.0 => NaN)
    Ygrids = [Ydistrib[:, col] for col in 1:size(Ydistrib,2)]
    Ygrids = reshape.(Ygrids, dims...)
    distributions = SimpleSDMResponse.(Ygrids, lims...)
    push!(layers, distributions)
end
distributions, prob_distrib, lower_distrib, upper_distrib = layers;
distributions

## Export results
# save_data = true
if (@isdefined save_data) && save_data == true
    @save joinpath("data", "jld2", "bart-distributions.jld2") distributions prob_distrib lower_distrib upper_distrib
    _zip_jld2(joinpath("data", "jld2", "bart-distributions.zip"),
              joinpath("data", "jld2", "bart-distributions.jld2"))
    touch(joinpath("data", "jld2", "bart-distributions.jld2"))
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])
