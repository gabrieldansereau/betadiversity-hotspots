#### Data preparation ####

message("Loading & preparing data")

# subset_qc <- TRUE

## 1. Load data ####

# Load full range data
spa_full <- read_tsv(here("data", "proc", "distributions_spa_full.csv"))
env_full <- read_tsv(here("data", "proc", "distributions_env_full.csv"))
spe    <- read_tsv(here("data", "proc", "distributions_spe_full.csv"))

# Load QC data (optional)
spa_qc <- read_tsv(here("data", "proc", "distributions_spa_qc.csv"))


## 2. Prepare data ####

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
  spe    <- dplyr::filter(spe, site %in% sites_qc)
  sites_obs <- dplyr::intersect(sites_qc, sites_obs)
  inds_obs  <- which(sites_qc %in% inds_obs) # not the same for QC data
}

# Filter spa & env to observed sites only
spa <- dplyr::filter(spa_full, site %in% sites_obs)
env <- dplyr::filter(env_full, site %in% sites_obs)
# Expand spe to full scale
spe_full <- as_tibble(
   matrix(
      NA_real_,
      nrow = nrow(spa_full), ncol = ncol(spe),
      dimnames = list(NULL, names(spe))
   )
)
spe_full[inds_obs,] <- spe

# Remove site with NAs for landcover variables
(inds_withNAs <- unique(unlist(map(env, ~ which(is.na(.x))))))
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
vars <- dplyr::select(vars, -site)

# Remove species without observations
(spe_withoutobs <- names(which(colSums(spe) == 0)))
if (length(spe_withoutobs > 0)) {
  spe <- dplyr::select(spe, -all_of(spe_withoutobs))
  spe_full <- dplyr::select(spe_full, -all_of(spe_withoutobs))
}
