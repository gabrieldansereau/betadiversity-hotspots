import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

## Conditional arguments
# save_figures = true

## Load data
#=
spe = CSV.read(joinpath("data", "proc", "distributions_spe_full.csv"), header=true, delim="\t")
spa = CSV.read(joinpath("data", "proc", "distributions_spa_full.csv"), header=true, delim="\t")
env = CSV.read(joinpath("data", "proc", "distributions_env_full.csv"), header=true, delim="\t")
=#

## Perform RandomForests
using RCall
# @rput spe spa env
begin
    R"""
    ## 0. Load packages ####
    library(embarcadero)
    library(tidyverse)
    library(viridis)
    library(furrr)


    ## 1. Load data ####

    # Load data
    spa_full <- read_tsv("data/proc/distributions_spa_full.csv")
    env_full <- read_tsv("data/proc/distributions_env_full.csv")

    # Subset to QC data (optional)
    spa_full <- read_tsv("data/proc/distributions_spa_qc.csv")
    sites_qc <- spa_full$site
    env_full <- filter(env_full, site %in% sites_qc)

    # Combine environmental variables
    vars_full <- left_join(spa_full, env_full, by = "site")

    # Select fewer variables
    xnames <- c("lat", "lon", "wc1", "wc12", paste0("lc", c(1:5, 7:10))) # lc6 always zero
    vars_df <- vars_full[,xnames]
 
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

    ## 2. Models ####
    
    # Load models
    load("data/proc/bart_models_qc.RData")

    # Quantile Predictions
    system.time(
        predictions <- future_map(
            sdms,
            function(x) predict(
                object = x, 
                x.layers = vars_stack,
                quantiles = c(0.025, 0.975),
                splitby = 20
            )
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
    summary_inner <- function(object){
        # Fit
        if(class(object)=='bart') {
            fitobj <- object$fit
        } else if(class(object)=='rbart') {
            fitobj <- object$fit[[1]] 
        }
        
        # Call
        (objcall <- object$call)
        
        # Predictors
        (predictors <- paste(attr(fitobj$data@x, "term.labels"), sep=' '))
        
        # Labels
        (true.vector <- fitobj$data@y)
        
        # Predictions
        (predictions <- colMeans(pnorm(object$yhat.train)))
        
        # Prediction instance
        (pred <- prediction(predictions, true.vector))
        # Performance instance
        (perf.tss <- performance(pred, "sens", "spec"))
        # TSS values list
        tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
        # TSS values ~ threshold dataframe
        (tss.df <- tibble(alpha=perf.tss@alpha.values[[1]],tss=tss.list))
        
        # AUC
        (auc <- performance(pred,"auc")@y.values[[1]])
        
        # Threshold
        (thresh <- min(tss.df$alpha[which(tss.df$tss==max(tss.df$tss))]))
        
        # TSS
        (tss <- tss.df[which(tss.df$alpha==thresh),'tss'][[1]])
        
        # Type I error rate (false positive)
        (type_I <-  1 - perf.tss@y.values[[1]][which(perf.tss@alpha.values[[1]] == thresh)])
        # Type II error rate (false negative)
        (type_II <- 1 - perf.tss@x.values[[1]][which(perf.tss@alpha.values[[1]] == thresh)])
        
        diagnostics <- list(
            fit = fitobj,
            call = objcall,
            predictors = predictors,
            labels = true.vector,
            predictions = predictions,
            auc = auc,
            threshold = thresh,
            tss = tss,
            type_I = type_I,
            type_II = type_II
        )
        return(diagnostics)
    }
    summaries <-  map(sdms, summary_inner)

    # Organize as tibble
    results <- tibble(
        spe = names(spe),
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
Y = replace(predictions_full, missing => NaN)
# Set values to NaN if no species present
inds_zeros = findall(map(x -> all(iszero.(x)), eachrow(Y)))
Y[inds_zeros,:] .= NaN
# Get Yobs
inds_obs = findall(map(x -> !any(isnan.(x)), eachrow(Y)))
inds_notobs = findall(map(x -> any(isnan.(x)), eachrow(Y)))
Yobs = Y[inds_obs,:]
# Apply Hellinger transformation
@rput Yobs
begin
    R"""
        library(vegan)
        Ytransf <- decostand(Yobs, "hel")
    """
end
@rget Ytransf

## Create distributions

# Load raw distributions (for grid size)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions

# Create RF distribution layers
Ydistrib = replace(Y, 0.0 => NaN)
rf_grids = [reshape(Ydistrib[:,i], size(distributions[1].grid)) for i in 1:size(Ydistrib, 2)]
# map(x -> reshape(Ydistrib[:,x], size(distributions[1]), 1:size(Ydistrib, 2)))
distributions = SimpleSDMResponse.(rf_grids, distributions[1].left, distributions[1].right,
                                         distributions[1].bottom, distributions[1].top)

## Export results
# save_data = true
if @isdefined save_data && save_data == true
    @save joinpath("data", "jld2", "rf-distributions.jld2") distributions
end

## Get richness & LCBD

richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])
