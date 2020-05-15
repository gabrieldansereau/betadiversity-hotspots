## 0. Load packages ####
library(tidyverse)
library(embarcadero)


## 1. Load data ####

# Load data
spe_full <- read.csv("data/proc/distributions_spe_full.csv", header=TRUE, sep="\t")  %>% as_tibble() 
spa_full <- read.csv("data/proc/distributions_spa_full.csv", header=TRUE, sep="\t")  %>% as_tibble()
env_full <- read.csv("data/proc/distributions_env_full.csv", header=TRUE, sep="\t")  %>% as_tibble()

# Select observed sites only
inds_obs <- spe_full$site

# Subset to QC data (optional)
spa_full <- read.csv("data/proc/distributions_spa_qc.csv", header=TRUE, sep="\t")  %>% as_tibble()
env_full <- env_full[spa_full$site,]
inds_obs <- intersect(spa_full$site, inds_obs)
spe_full <- spe_full[spe_full$site %in% inds_obs,]

# Filter datasets to observed sites only
spa <- subset(spa_full, site %in% inds_obs)
env <- subset(env_full, site %in% inds_obs)
spe <- spe_full

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
spe <- subset(spe, select = -site)
env <- subset(env, select = -site)
spa <- subset(spa, select = -site)
vars <- subset(vars, select = -site)

# Remove species without observations
(inds_withoutobs <- c(which(sapply(spe, sum) == 0)))
if (length(inds_withoutobs > 0)) {
    spe <- subset(spe, select = -inds_withoutobs)
}

# Select species with ~ same number of presences & absences
colSums(spe)/nrow(spe) # 17 seems good
sp <- subset(spe, select = "sp17") # black-throated blue warbler
table(sp)

# Select fewer variables
xnames <- c("lat", "lon", "wc1", "wc12", paste0("lc", 1:10))

## 2. Basic BART model ####

# Train model
system.time(
    sdm <- bart(y.train = sp$sp17,
                x.train = vars[,xnames],
                keeptrees = TRUE
                )
)

# Model diagnostics
summary(sdm)

# Create raster layer
df_to_layer <- function(x, lons, lats){
    mat <- matrix(data = x, nrow = uniqueN(lats), ncol = uniqueN(lons))
    layer <- raster(mat[nrow(mat):1,],
        xmn=min(lons), xmx=max(lons), 
        ymn=min(lats), ymx=max(lats)
        )
    return(layer)
}
wc_layer <- df_to_layer(x = vars_full$wc1, lons = vars_full$lon, lats = vars_full$lat)
plot(wc_layer)

# Stack raster layers
vars_df <- as.data.frame(vars_full[,xnames])
vars_layers <- lapply(1:ncol(vars_full[,xnames]), function(x) df_to_layer(vars_df[,x], lons = vars_df$lon, lats = vars_df$lat)
vars_stack <- stack(vars_layers)
names(vars_stack) <- xnames
vars_stack

# Predict species distribution
map <- predict(object = sdm, x.layers = vars_stack, splitby = 20, quiet = TRUE)

# Plot predictions
plot(map, main='Predicted probability', 
     box=F, axes=F)
