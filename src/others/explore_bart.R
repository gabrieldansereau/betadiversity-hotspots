## 0. Load packages ####
library(embarcadero)
library(tidyverse)
library(viridis)
library(furrr)
plan(multiprocess)


## 1. Load data ####

# Load data
spa_full <- read_tsv("data/proc/distributions_spa_full.csv")
env_full <- read_tsv("data/proc/distributions_env_full.csv")
spe      <- read_tsv("data/proc/distributions_spe_full.csv") 

# Select observed sites only
sites_obs <- spe$site
inds_obs  <- which(spa_full$site %in% sites_obs)
all(sites_obs == inds_obs) # happen to be the same at full scale

# Subset to QC data (optional)
spa_full <- read_tsv("data/proc/distributions_spa_qc.csv")
sites_qc <- spa_full$site
env_full <- filter(env_full, site %in% sites_qc)
spe      <- filter(spe, site %in% sites_qc)
sites_obs <- intersect(sites_qc, sites_obs)
inds_obs  <- which(sites_qc %in% inds_obs) # not the same for QC data

# Filter spa & env to observed sites only
spa <- filter(spa_full, site %in% sites_obs)
env <- filter(env_full, site %in% sites_obs)
# Expand spe to full scale
spe_full <- as_tibble(
     matrix(
          NaN, 
          nrow = nrow(spa_full), ncol = ncol(spe), 
          dimnames = list(NULL, names(spe))
     )
)
spe_full[inds_obs,] <- spe

# Remove site with NAs for landcover variables
(inds_withNAs <- unique(unlist(sapply(env, function(x) which(is.na(x))))))
if (length(inds_withNAs) > 0) {
    spe <- spe[-inds_withNAs,]
    spa <- spa[-inds_withNAs,]
    env <- env[-inds_withNAs,]
}

# Combine environmental variables
vars <- left_join(spa, env, by = "site")
vars_full <- left_join(spa_full, env_full, by = "site")
# Remove site column
spe  <- select(spe, -site)
env  <- select(env, -site)
spa  <- select(spa, -site)
vars <- select(vars, -site)

# Remove species without observations
(inds_withoutobs <- c(which(sapply(spe, sum) == 0)))
if (length(inds_withoutobs > 0)) {
    spe <- subset(spe, select = -inds_withoutobs)
    spe_full <- subset(spe_full, select = -inds_withoutobs)
}

# Select 1 species with ~ same number of presences & absences
colSums(spe)/nrow(spe) # 17 seems good
sp <- select(spe, sp17) # black-throated blue warbler
sp_full <- select(spe_full, sp17) # black-throated blue warbler
table(sp)

# Select fewer variables
xnames <- c(paste0("wc", c(1:19)), paste0("lc", c(1:5, 7:10))) # lc6 always zero

## 2. Create layers ####

# Create raster layer
df_to_layer <- function(x, lons, lats){
    mat <- matrix(data = x, nrow = uniqueN(lats), ncol = uniqueN(lons))
    layer <- raster(
        mat[nrow(mat):1,],
        xmn=min(lons), xmx=max(lons), 
        ymn=min(lats), ymx=max(lats)
     )
    return(layer)
}
wc_layer <- df_to_layer(x = vars_full$wc1, lons = vars_full$lon, lats = vars_full$lat)
sp_layer <- df_to_layer(x = sp_full[[1]], lons = vars_full$lon, lats = vars_full$lat)
plot(wc_layer, col = "grey", legend = FALSE)
plot(sp_layer, main = "Setophaga caerulescens", col = viridis(2), add = TRUE)

# Stack raster layers
vars_df <- vars_full[,xnames]
vars_layers <- map(
    vars_df, 
    function(x) df_to_layer(x, lons = vars_full$lon, lats = vars_full$lat)
)
(vars_stack <- stack(vars_layers, names = xnames))
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
# Extract inner values
source("src/lib/R/bart.R")
diagnostics <- summary_inner(sdm)

# Predict species distribution
predictions <- predict(object = sdm, x.layers = vars_stack, splitby = 20)

# Plot predictions
par(mfrow=c(1,2), mar=c(10,4,10,7))
plot(wc_layer, main = "Setophaga caerulescens", col = "grey", 
    legend=FALSE, box = FALSE, axes = F)
plot(sp_layer, add = TRUE, main = "Setophaga caerulescens", col = viridis(2))
plot(predictions[[1]], main='Predicted probability', col = viridis(255),
     box=F, axes=F)
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)

# Thresholded predictions
threshold <- diagnostics$threshold
par(mfrow=c(1,2), mar=c(10,4,10,7))
plot(wc_layer, main = "Setophaga caerulescens", col = "grey", 
     legend=FALSE, box = FALSE, axes = F)
plot(sp_layer, add = TRUE, main = "Setophaga caerulescens", col = viridis(2))
plot(predictions[[1]] > threshold, main='Predicted outcome', col = viridis(255),
     box=F, axes=F)
par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)


## 3. More advanced model ####

# Lower-upper quantiles
pred_quant <- predict(sdm, vars_stack, quantiles=c(0.025, 0.975), splitby = 20)

# Plot quantiles
par(mfrow=c(2,2))
par(mar=c(2,1,3,7))
plot(pred_quant[[1]], 
     main = 'Posterior mean', 
     box=F, axes=F, col = viridis(255))
plot(pred_quant[[2]], 
     main = 'Lower 95% CI bound', 
     box=F, axes=F, col = viridis(255))
plot(pred_quant[[3]], 
     main = 'Upper 95% CI bound', 
     box=F, axes=F, col = viridis(255))
plot(pred_quant[[3]]-pred_quant[[2]], 
     main = 'Credible interval width', 
     box=F, axes=F, col = viridis(255))

# Plot quantiles with threshold
plot(pred_quant[[1]] > threshold, 
     main = 'Posterior mean', 
     box=F, axes=F, col = viridis(255))
plot(pred_quant[[2]] > threshold, 
     main = 'Lower 95% CI bound', 
     box=F, axes=F, col = viridis(255))
plot(pred_quant[[3]] > threshold, main = 'Upper 95% CI bound', 
     box=F, axes=F, col = viridis(255))
plot((pred_quant[[3]]> threshold) - (pred_quant[[2]] > threshold),
     main = 'Credible interval width', 
     box=F, axes=F, col = viridis(255)) # Thresholded presence difference
plot(pred_quant[[3]] - pred_quant[[2]] > threshold,
     main = 'Credible interval width', 
     box=F, axes=F, col = viridis(255)) # Probability difference vs threshold
par(mfrow=c(1,1))

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

## 5. Variable selection ####

varimp(sdm, plot = TRUE)
system.time(
    varimp.diag(vars[,xnames], sp[[1]], iter = 50)
) # ~ 15 min

# Stepwise variable set reduction
system.time(
    step.model <- variable.step(
        x.data = vars[,xnames], 
        y.data = sp[[1]]
    )
) # ~5 min
step.model
# "wc1"  "wc2"  "wc7"  "wc8"  "wc10" "wc11" "lc2"  "lc3"  "lc5" 

## 6. Partial dependence plots

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

## 7. Multi-species models ####

# Run for all species
set.seed(42)
system.time(
    sdms <- future_map(spe,
        function(x) bart(
            y.train = x,
            x.train = vars[,xnames],
            keeptrees = TRUE
        )
    )
) # ~ 4 min., 45 sec. in parallel

# Extract summary statistics
summary(sdms[[1]])
summaries <-  map(sdms, summary_inner)
summaries[[20]]
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
summary(results)

# Export results
# save(sdms, file = "data/proc/bart_models_qc.RData")
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

# Plot predictions
pred_plot <- ggplot(spa_full, aes(lon, lat)) +
    geom_raster(aes(fill = pred_df$sp17)) +
    scale_fill_viridis_c(na.value = "white") +
    ggtitle("Predictions") + 
    coord_quickmap() +
    theme_minimal()
pred_plot
pred_plot + geom_raster(aes(fill = lower_df$sp17))
pred_plot + geom_raster(aes(fill = upper_df$sp17))
pred_plot + geom_raster(aes(fill = upper_df$sp17 - lower_df$sp17))
pred_plot + geom_raster(aes(fill = pres_df$sp17))
