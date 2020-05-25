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

# Check prepared data
spa_full
env_full
spe_full
vars_full
spa
env
spe

# Select 1 species with ~ same number of presences & absences
colSums(spe)/nrow(spe) # 17 seems good
sp <- select(spe, sp17) # black-throated blue warbler
sp_full <- select(spe_full, sp17) # black-throated blue warbler
table(sp)

# Select fewer variables
# xnames <- c("wc1", "wc12", paste0("lc", c(1:5, 7:10))) # lc6 always zero
# xnames <- c(paste0("wc", c(1:19)), paste0("lc", c(1:5, 7:10))) # lc6 always zero
# xnames <- c("wc1", "wc2", "wc7", "wc8", "wc10", "wc11", "lc2", "lc3", "lc5") # (all vars - QC)
# xnames <- c("lat", "lon", "wc1", "wc2", "wc8",  "wc11", "lc2", "lc5", "lc8") # (all vars - QC with lat-lon)
# xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:5, 7:10))) # ~ vars from CCHF vignette + landcover
# xnames <- c("wc1", "wc2", "wc5", "wc6", "wc12", "wc13", "wc14", "wc15", "lc2", "lc3", "lc5", "lc8") # stepwise selection on CCHF vignette
# xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:10)))
xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:3,5,7:10)))

## 2. Create layers ####

# Create raster layers
sp_layer <- df_to_layer(x = sp_full[[1]], lons = vars_full$lon, lats = vars_full$lat)
vars_layers <- map(
    vars_full[,xnames], 
    function(x) df_to_layer(x, lons = vars_full$lon, lats = vars_full$lat)
)
wc_layer <- vars_layers$wc1

# Stack variables layers
vars_stack <- stack(vars_layers, names = xnames)
plot(vars_stack)


## 3. Basic BART model ####

# Train model
set.seed(42)
system.time(
    sdm <- bart(
        y.train = sp[[1]],
        x.train = vars[,xnames],
        keeptrees = TRUE
    )
) # 5 sec. for QC, 90 sec. at full scale

# Model diagnostics
summary(sdm)
varimp(sdm, plot = TRUE)
# Extract inner values
diagnostics <- summary_inner(sdm)

# Predict species distribution
predictions <- predict2.bart(sdm, vars_stack, quantiles=c(0.025, 0.975), splitby = 20)

# Plot probability predictions
plot(
    wc_layer, 
    main = "Original distribution",
    col = "grey",  
    legend = FALSE, box = FALSE, axes = FALSE
)
plot(sp_layer, add = TRUE, col = viridis(2))
plot(
    predictions[[1]], 
    main = 'Probability predictions - Posterior mean', 
    col = viridis(255),
    zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)
plot(
    predictions[[2]], 
    main = 'Probability predictions - Lower 95% CI bound', 
    col = viridis(255),
    zlim = c(0, 1), 
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)
plot(
    predictions[[3]], 
    main = 'Probability predictions - Upper 95% CI bound', 
    col = viridis(255),
    zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)
plot(
    predictions[[3]] - predictions[[2]], 
    main = 'Probability predictions - Credible interval width', 
    col = viridis(255),
    zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)

# Plot threshold predictions
threshold <- diagnostics$threshold
plot(
    predictions[[1]] > threshold, 
    main = 'Threshold predictions - Posterior mean', 
    col = viridis(2),
    box = FALSE, axes = FALSE
)
plot(
    predictions[[2]] > threshold, 
    main = 'Threshold predictions - Lower 95% CI bound', 
    col = viridis(2),
    box = FALSE, axes = FALSE
)
plot(
    predictions[[3]] > threshold, 
    main = 'Threshold predictions - Upper 95% CI bound', 
    col = viridis(255),
    box = FALSE, axes = FALSE
)
plot(
    (predictions[[3]] > threshold) - (predictions[[2]] > threshold), 
    main = 'Threshold predictions - Upper-lower site difference', 
    col = viridis(2),
    box = FALSE, axes = FALSE
)


## 4. Variable selection ####

varimp(sdm, plot = TRUE)
system.time(
    # varimp.diag(vars[,xnames], sp[[1]], iter = 50)
) # ~ 15 min

# Stepwise variable set reduction
system.time(
    step.model <- variable.step(
        x.data = vars[,xnames], 
        y.data = sp[[1]]
    )
) # ~5 min
step.model

# 2:30:00 for full scale 
# wc1  wc5  wc6  wc12 wc14 wc15 lc7  lc8  lc9

## 5. Partial dependence plots ####

system.time(
    partial(
        sdm, 
        x.vars = c('wc1'),
        trace = FALSE,
        ci = FALSE,
        smooth = 5,
        equal = TRUE
    )
) # ~ 3 min.

# Spartial dependence plot
system.time(
    p <- spartial(sdm, vars_stack, x.vars='wc1', equal=TRUE, smooth=5)
) # ~ 3 min.
plot(p)


## 6. Multi-species model training ####

# Run for all species
set.seed(42)
system.time(
    sdms <- pbapply::pblapply(
        spe[1:10],
        function(x) bart(
            y.train = x,
            x.train = vars[,xnames],
            keeptrees = TRUE,
            verbose = FALSE
        )
    )
) # ~ 4 min., 45 sec. in parallel

# Export results
# save_models <- TRUE
if (exists("save_models") && isTRUE(save_models)) {
    message("Saving models to RData")
    save(sdms, file = "data/proc/bart_models.RData")
} else {
    message("Loading models from RData")
    load("data/proc/bart_models.RData")
}

# Extract summary statistics
summary(sdms[[1]])
summaries <-  map(sdms, summary_inner)
summaries[[2]]
str(summaries)

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


## 7. Multi-species predictions ####

# Quantile Predictions
sdms <- sdms[1:10] 
system.time(
    predictions <- future_map(
        sdms[1:2],
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

# Plot predictions
pred_plot <- ggplot(spa_full, aes(lon, lat)) +
    geom_raster(aes(fill = pred_df$sp1)) +
    scale_fill_viridis_c(na.value = "white") +
    ggtitle("Predictions") + 
    coord_quickmap() +
    theme_minimal()
pred_plot
pred_plot + geom_raster(aes(fill = lower_df$sp1))
pred_plot + geom_raster(aes(fill = upper_df$sp1))
pred_plot + geom_raster(aes(fill = upper_df$sp1 - lower_df$sp1))
pred_plot + geom_raster(aes(fill = pres_df$sp1))


## 8. Others ####

# Show first iterations
system.time(
     plot.mcmc(sdm, vars_stack, iter=5, quiet = TRUE)
)

# Show burn-in with fewer trees
set.seed(42)
sdm.tiny <- bart(
    y.train=sp[[1]],
    x.train=vars[,xnames],
    keeptrees = TRUE,
    ntree=5, # 5 tree models
    nskip=0
) # No burnin
summary(sdm.tiny)
plot.mcmc(sdm.tiny, vars_stack, iter=100)

# Timelapse of tree learning
library(animation)
saveGIF(
    plot.mcmc(sdm, vars_stack, iter=50), 
    movie.name = "Timelapse.gif", #interval = 0.15, 
    ani.width = 800, ani.height = 400
)

