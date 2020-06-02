library(conflicted)
library(tidyverse)
library(here)
library(embarcadero)
library(viridis)
library(furrr)
plan(multiprocess)
# Resolve conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("select", "dplyr")
# Custom functions
source(here("src", "lib", "R", "bart.R"))

## 1. Load data ####

message("Loading & preparing data")

# Load data
spa_full <- read_tsv(here("data", "proc", "distributions_spa_full.csv"))
env_full <- read_tsv(here("data", "proc", "distributions_env_full.csv"))
spe      <- read_tsv(here("data", "proc", "distributions_spe_full.csv"))

# Load QC data (optional)
spa_qc <- read_tsv(here("data", "proc", "distributions_spa_qc.csv"))

# Prepare data
# subset_qc <- TRUE # subset to QC data (optional)
source(here("src", "02_training_data-preparation.R"))

# Remove site with NAs for landcover variables
if (length(inds_withNAs) > 0) {
    richness_values <- richness_values[-inds_withNAs]
    lcbd_values <- lcbd_values[-inds_withNAs]
}

# Select variables
xnames <- select(vars, -c(lat, lon)) %>% names()

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


## 3. Model training ####

# Group richness & lcbd
values_df <- tibble(
    richness = richness_values,
    lcbd = lcbd_values
)

# Function to run BART in parallel
bart_parallel <- function(...) {
    # BART model
    model <- bart(...)
    # Touch state so that saving will work
    invisible(model$fit$state)
    return(model)
}

# Train model
set.seed(42)
system.time(
    models <- future_map(
        values_df,
        ~ bart_parallel(
            y.train = .x,
            x.train = vars[,xnames],
            keeptrees = TRUE
        )
    )
) # 90 sec.
varimp(models[[1]], plot = TRUE)
varimps <-  map_dfr(models, varimp, .id = "value") %>% 
    rename(vars = names)
varimps %>% 
    pivot_wider(names_from = "value", values_from = "varimps") %>% 
    mutate(diff = richness - lcbd)
ggplot(varimps, aes(vars)) + 
    geom_bar(aes(weight = varimps, fill = value)) #, position = "dodge2")

ggplot(varimps, aes(vars, varimps)) + 
    geom_boxplot(aes(colour = value))

## 4. Predictions ####
system.time(
    predictions <- future_map(
        models,
        ~ predict2.bart(
            .x,
            vars_stack,
            quantiles = c(0.025, 0.0975),
            splitby = 20
        )
    )
) # 2 min.

# Plot richness
plot(
    predictions$richness[[1]], 
    main = 'Probability predictions - Posterior mean', 
    col = viridis(255),
    # zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)
# Negative richness??
plot(
    predictions$richness[[1]] < 1, 
    main = 'Probability predictions - Posterior mean', 
    col = viridis(255),
    # zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)
plot(
    predictions$lcbd[[1]], 
    main = 'Probability predictions - Posterior mean', 
    col = viridis(255),
    zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)

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