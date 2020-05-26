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
# save_models <- TRUE # save & overwrite models


## 1. Load data ####

# Load data
spa_full <- read_tsv("data/proc/distributions_spa_full.csv")
env_full <- read_tsv("data/proc/distributions_env_full.csv")
spe      <- read_tsv("data/proc/distributions_spe_full.csv") 

# Load QC data (optional)
spa_qc <- read_tsv("data/proc/distributions_spa_qc.csv")

# Prepare data
# subset_qc <- TRUE # subset to QC data (optional)
source("src/02_training_data-preparation.R")

# Select 1 species with ~ same number of presences & absences
colSums(spe)/nrow(spe) # 17 seems good
sp <- select(spe, sp17) # black-throated blue warbler
sp_full <- select(spe_full, sp17) # black-throated blue warbler
table(sp)

# Select fewer variables
xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:3,5,7:10)))

## 2. Create layers ####

# Create raster layers
vars_layers <- map(
    vars_full[,xnames], 
    function(x) df_to_layer(x, lons = vars_full$lon, lats = vars_full$lat)
)
wc_layer <- vars_layers$wc1

# Stack variables layers
vars_stack <- stack(vars_layers, names = xnames)
plot(vars_stack)


## 3. Multi-species model training ####

bart_parallel <- function(x.train, y.train, ...) {
    sdm <- bart(
        y.train = y.train,
        x.train = x.train,
        keeptrees = TRUE,
        verbose = FALSE,
        ...
    )
    invisible(sdm$fit$state)
    return(sdm)
}

# Run for all species
# create_models <- TRUE
if (exists("create_models") && isTRUE(save_models)){
    set.seed(42)
    system.time(
        sdms <- future_map(
            spe[1:2],
            function(x) bart_parallel(
                y.train = x,
                x.train = vars[,xnames]
            )
        )
    ) # ~ 4 min., 45 sec. in parallel
    invisible(map(sdms, function(x) x$fit$state))
    sdms_backup <- sdms
}

# Export results
# save_models <- TRUE
if (exists("save_models") && isTRUE(save_models)) {
    message("Saving models to RData")
    save(sdms, file = "data/proc/bart_models_01-02.RData")
} else {
    message("Loading models from RData")
    load("data/proc/bart_models_01-02.RData")
}

# Extract summary statistics
summaries <-  future_map(sdms, summary_inner)
summaries[[2]]

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
print(results, n = Inf)
summary(results)

# Extract variable importance
varimps <- map(sdms, varimp) %>% 
    map_df(~ .$varimps) %>% 
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
) # ~ 8 min., 1 min. in parallel
predictions_orig <- predictions
predictions_touch <- predictions

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

# Presence-absence dataframe
pres_df <- map2_df(
    pred_df, results$threshold, 
    function(pred, thresh) ifelse(pred > thresh, 1, 0) 
)
pres_df

# Plot predictions
pred_plot <- ggplot(spa_full, aes(lon, lat)) +
    scale_fill_viridis_c(na.value = "white", "Value") +
    coord_quickmap() +
    theme_minimal()
sp_no <- 2
pred_plot + geom_raster(aes(fill = pred_df[[sp_no]])) + ggtitle("Posterior mean")
pred_plot + geom_raster(aes(fill = lower_df[[sp_no]])) + ggtitle("Lower CI")
pred_plot + geom_raster(aes(fill = upper_df[[sp_no]])) + ggtitle("Upper CI")
pred_plot + geom_raster(aes(fill = upper_df[[sp_no]] - lower_df[[sp_no]])) + ggtitle("CI difference")
pred_plot + geom_raster(aes(fill = pres_df[[sp_no]])) + ggtitle("Presence-absence")
