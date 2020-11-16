## 0. Load packages ####
source(file.path("src", "required.R"))

# Conditional evaluations
subset_qc <- TRUE # subset to QC data (optional)
# save_models <- TRUE # save & overwrite models


## 1. Load data ####

# Load data
source(here("src", "02a_training_data-preparation.R"))

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
    ~ df_to_layer(.x, lons = vars_full$lon, lats = vars_full$lat)
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


## 6. Others ####

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

