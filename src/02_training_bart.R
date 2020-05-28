## 0. Load packages ####
library(conflicted)
library(embarcadero)
library(tidyverse)
library(viridis)
library(furrr)
plan(multiprocess)

# Resolve conflicts
conflict_scout()
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("select", "dplyr")

# Custom functions
source("src/lib/R/bart.R")

# Conditional evaluations
# subset_qc <- TRUE # subset to QC data (optional)
# create_models <- TRUE # train models
# save_models <- TRUE # save & overwrite models


## 1. Load data ####

message("Loading & preparing data")

# Load data
spa_full <- read_tsv("data/proc/distributions_spa_full.csv")
env_full <- read_tsv("data/proc/distributions_env_full.csv")
spe      <- read_tsv("data/proc/distributions_spe_full.csv") 

# Load QC data (optional)
spa_qc <- read_tsv("data/proc/distributions_spa_qc.csv")

# Prepare data
# subset_qc <- TRUE # subset to QC data (optional)
source("src/02_training_data-preparation.R")

# Select fewer variables
xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:3,5,7:10)))

## 2. Create layers ####

message("Creating layers")

# Create raster layers
vars_layers <- map(
    vars_full[,xnames], 
    ~ df_to_layer(.x, lons = vars_full$lon, lats = vars_full$lat)
)
wc_layer <- vars_layers$wc1

# Stack variables layers
vars_stack <- stack(vars_layers, names = xnames)
plot(vars_stack)


## 3. Multi-species model training ####

message("Training multi-species models")

# Split species in groups
# max_procs <- availableCores() - 1
max_procs <- 7
spe_num <- 1:ncol(spe)
spe_splits <- split(spe_num, ceiling(spe_num/max_procs))
spe_groups <- map(spe_splits, ~ select(spe, all_of(.)))

# Function to run BART in parallel
bart_parallel <- function(x.train, y.train, ...) {
    # BART SDM
    sdm <- bart(
        y.train = y.train,
        x.train = x.train,
        keeptrees = TRUE,
        verbose = FALSE,
        ...
    )
    # Touch state so that saving will work
    invisible(sdm$fit$state)
    return(sdm)
}

# Prepare global sets
sdms_list <- list()
predictions_list <- list()
results_global <- tibble(
    spe = names(spe),
    auc = NaN,
    threshold = NaN,
    tss = NaN,
    type_I = NaN,
    type_II = NaN
)
varimps_global <- matrix(
    NaN,
    nrow = length(xnames), ncol = ncol(spe),
    dimnames = list(NULL, names(spe))
) %>% 
    as_tibble() %>% 
    mutate(vars = xnames) %>% 
    select(vars, everything())
pred_df_global <- matrix(
    NaN,
    nrow = nrow(vars_full), ncol = ncol(spe),
    dimnames = list(NULL, names(spe))
) %>% as_tibble()
lower_df_global <- pred_df_global
upper_df_global <- pred_df_global
pres_df_global <- pred_df_global


# Run for each group
system.time(
for (gp in seq_along(spe_groups)) {

    message(paste0("Training multi-species group (", gp, "/", length(spe_groups)), ")")

    # create_models <- TRUE
    if (exists("create_models") && isTRUE(create_models)){
        set.seed(42)
        system.time(
            sdms <- future_map(
                spe_groups[[gp]],
                ~ bart_parallel(
                    y.train = .x,
                    x.train = vars[,xnames]
                ),
                .progress = TRUE
            )
        ) # ~ 4 min in parallel
    }

    # Export results
    # save_models <- TRUE
    filepath <- paste0("data/proc/bart_models", gp, ".RData")
    if (exists("save_models") && isTRUE(save_models)) {
        message("Saving models to RData")
        save(sdms, file = filepath)
    } else {
        message("Loading models from RData")
        load(filepath)
    }

    # Extract summary statistics
    summaries <-  future_map(sdms, summary_inner)

    # Organize as tibble
    results <- tibble(
        spe = names(summaries),
        auc = map_dbl(summaries, "auc"),
        threshold = map_dbl(summaries, "threshold"),
        tss = map_dbl(summaries, "tss"),
        type_I = map_dbl(summaries, "type_I"),
        type_II = map_dbl(summaries, "type_II")
    )
    print(results, n = Inf)
    summary(results)

    # Extract variable importance
    varimps <- map(sdms, varimp) %>% 
        map_df("varimps") %>% 
        mutate(vars = xnames) %>% 
        select(vars, everything())
    varimps
    varimps %>% 
        pivot_longer(-vars, names_to = "spe", values_to = "varimp") %>% 
        pivot_wider(spe, names_from = "vars", values_from = "varimp")
    varimps %>% 
        pivot_longer(-vars, names_to = "spe", values_to = "varimp") %>% 
        group_by(vars) %>% 
        summarize(mean = mean(varimp)) %>% 
        arrange(desc(mean))


    ## 4. Multi-species predictions ####

    message("Predicting species distributions")

    # Quantile Predictions
    system.time(
        predictions <- future_map(
            sdms,
            # sdms_backup,
            function(x) predict2.bart(
                object = x, 
                x.layers = vars_stack,
                quantiles = c(0.025, 0.975),
                splitby = 20,
                quiet = TRUE
            ),
            .progress = TRUE
        )
    ) # ~ 3 min in parallel

    # Predictions
    pred_df <- predictions %>% 
        map(~ .x$layer.1) %>% 
        stack() %>% 
        as.data.frame(xy = TRUE) %>% 
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    pred_df
    # Lower quantiles
    lower_df <- predictions %>% 
        map(~ .x$layer.2) %>% 
        stack() %>% 
        as.data.frame(xy = TRUE) %>% 
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    lower_df
    # Upper quantiles
    upper_df <- predictions %>% 
        map(~ .x$layer.3) %>%
        stack() %>% 
        as.data.frame(xy = TRUE) %>%  
        as_tibble() %>% 
        arrange(x, y) %>% 
        select(-c(x, y))
    upper_df

    # Presence-absence dataframe
    pres_df <- map2_df(
        pred_df, results$threshold, 
        function(pred, thresh) ifelse(pred > thresh, 1, 0) 
    )
    pres_df

    # Export group results
    # sdms_list[[gp]] <- sdms
    # predictions_list[[gp]] <- predictions
    results_global[spe_splits[[gp]],] <- results
    varimps_global[spe_splits[[gp]]+1] <- varimps[,-1]
    pred_df_global[spe_splits[[gp]]] <- pred_df
    lower_df_global[spe_splits[[gp]]] <- lower_df
    upper_df_global[spe_splits[[gp]]] <- upper_df
    pres_df_global[spe_splits[[gp]]] <- pres_df
}
)

# Export to CSV
# save_predictions <- TRUE
if (exists("save_predictions") && isTRUE(save_predictions)) {
    write_tsv(results_global, "data/proc/bart_summaries.csv")
    write_tsv(varimps_global, "data/proc/bart_varimps.csv")
    write_tsv(pred_df_global, "data/proc/bart_predictions_prob.csv")
    write_tsv(lower_df_global, "data/proc/bart_predictions_lower.csv")
    write_tsv(upper_df_global, "data/proc/bart_predictions_upper.csv")
    write_tsv(pres_df_global, "data/proc/bart_predictions_pres.csv")
}

## 5. Visualize results ####

# Empty canvas
pred_plot <- ggplot(spa_full, aes(lon, lat)) +
    scale_fill_viridis_c(na.value = "white", "Value") +
    coord_quickmap() +
    theme_minimal()

# Plot predictions
sp_no <- 2
pred_plot + geom_raster(aes(fill = pred_df_global[[sp_no]])) + ggtitle("Posterior mean")
pred_plot + geom_raster(aes(fill = lower_df_global[[sp_no]])) + ggtitle("Lower CI")
pred_plot + geom_raster(aes(fill = upper_df_global[[sp_no]])) + ggtitle("Upper CI")
pred_plot + geom_raster(aes(fill = upper_df_global[[sp_no]] - lower_df_global[[sp_no]])) + ggtitle("CI difference")
pred_plot + geom_raster(aes(fill = pres_df_global[[sp_no]])) + ggtitle("Presence-absence")
