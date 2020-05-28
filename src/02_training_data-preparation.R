## Data preparation ####

# Select observed sites only
sites_obs <- spe$site
inds_obs  <- which(spa_full$site %in% sites_obs)
all(sites_obs == inds_obs) # happen to be the same at full scale

# Subset to QC data (optional)
# subset_qc <- TRUE
if (exists("subset_qc") && isTRUE(subset_qc)) {
    message("Subsetting to QC data")
    spa_full <- spa_qc
    sites_qc <- spa_full$site
    env_full <- dplyr::filter(env_full, site %in% sites_qc)
    spe      <- dplyr::filter(spe, site %in% sites_qc)
    sites_obs <- dplyr::intersect(sites_qc, sites_obs)
    inds_obs  <- which(sites_qc %in% inds_obs) # not the same for QC data
}

# Filter spa & env to observed sites only
spa <- dplyr::filter(spa_full, site %in% sites_obs)
env <- dplyr::filter(env_full, site %in% sites_obs)
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
spe  <- dplyr::select(spe, -site)
env  <- dplyr::select(env, -site)
spa  <- dplyr::select(spa, -site)
vars <- dplyr::select(vars, -site)

# Remove species without observations
(inds_withoutobs <- c(which(sapply(spe, sum) == 0)))
if (length(inds_withoutobs > 0)) {
    spe <- subset(spe, select = -inds_withoutobs)
    spe_full <- subset(spe_full, select = -inds_withoutobs)
}