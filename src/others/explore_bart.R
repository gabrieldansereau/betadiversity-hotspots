## 0. Load packages ####
library(embarcadero)
library(tidyverse)
library(viridis)


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
xnames <- c("lat", "lon", "wc1", "wc12", paste0("lc", c(1:5, 7:10))) # lc6 always zero

## 2. Basic BART model ####

# Train model
set.seed(42)
system.time(
    sdm <- bart(
        y.train = sp$sp17,
        x.train = vars[,xnames],
        keeptrees = TRUE
    )
)

# Model diagnostics
summary(sdm)

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
plot(wc_layer, col = "grey")
plot(sp_layer, main = "Setophaga caerulescens", col = viridis(2), add = TRUE)

# Stack raster layers
vars_df <- vars_full[,xnames]
vars_layers <- lapply(
    vars_df, 
    function(x) df_to_layer(x, lons = vars_df$lon, lats = vars_df$lat)
)
(vars_stack <- stack(vars_layers, names = xnames))
plot(vars_stack)
plot(vars_stack$wc1, main = "Temperature", col = "grey", legend = FALSE)
plot(sp_layer, main = "Setophaga caerulescens", col = viridis(2), add = TRUE)

# Predict species distribution
predictions <- predict(object = sdm, x.layers = vars_stack, splitby = 20, quiet = TRUE)

# Plot predictions
par(mfrow=c(1,2), mar=c(10,4,10,7))
plot(wc_layer, main = "Setophaga caerulescens", col = "grey", 
    legend=FALSE, box = FALSE, axes = F)
plot(sp_layer, add = TRUE, main = "Setophaga caerulescens", col = viridis(2))
plot(predictions[[1]], main='Predicted probability', col = viridis(255),
     box=F, axes=F)
par(mfrow=c(1,1))

# Thresholded predictions
summary(sdm) # Cutoff =  0.5556909
threshold <- 0.533461
par(mfrow=c(1,2), mar=c(10,4,10,7))
plot(wc_layer, main = "Setophaga caerulescens", col = "grey", 
     legend=FALSE, box = FALSE, axes = F)
plot(sp_layer, add = TRUE, main = "Setophaga caerulescens", col = viridis(2))
plot(predictions[[1]] > threshold, main='Predicted outcome', col = viridis(255),
     box=F, axes=F)
par(mfrow=c(1,1))

## 3. More advanced model ####

pred_quant <- predict(sdm, vars_stack, quantiles=c(0.025, 0.975), splitby = 20)

par(mfrow=c(2,2))
par(mar=c(2,1,3,7))
plot(pred_quant[[1]], 
     main = 'Posterior mean', 
     box=F, axes=F)
plot(pred_quant[[2]], 
     main = 'Lower 95% CI bound', 
     box=F, axes=F)
plot(pred_quant[[3]], 
     main = 'Upper 95% CI bound', 
     box=F, axes=F)
plot(pred_quant[[3]]-pred_quant[[2]], 
     main = 'Credible interval width', 
     box=F, axes=F)

plot(pred_quant[[1]] > threshold, 
     main = 'Posterior mean', 
     box=F, axes=F)
plot(pred_quant[[2]] > threshold, 
     main = 'Lower 95% CI bound', 
     box=F, axes=F)
plot(pred_quant[[3]] > threshold, main = 'Upper 95% CI bound', 
     box=F, axes=F)
plot((pred_quant[[3]]> threshold) - (pred_quant[[2]] > threshold),
     main = 'Credible interval width', 
     box=F, axes=F)
plot(pred_quant[[3]] - pred_quant[[2]] > threshold,
     main = 'Credible interval width', 
     box=F, axes=F)
par(mfrow=c(1,1))

# Show first iterations
plot.mcmc(sdm, vars_stack, iter=5, quiet = TRUE)

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
varimp.diag(vars[,xnames], sp[[1]], iter = 50)

# Stepwise variable set reduction
step.model <- variable.step(
    x.data = vars[,xnames], 
    y.data = sp[[1]]
)
step.model

## 6. Partial dependence plots

partial(
    sdm, 
    x.vars = c('wc1'),
    trace = FALSE,
    ci = FALSE,
    smooth = 5,
    equal = TRUE
)

# Spartial dependence plot
p <- spartial(sdm, vars_stack, x.vars='wc1', equal=TRUE, smooth=5)
plot(p)