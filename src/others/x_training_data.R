#### Prepare training data ####
# Exactly the same code as the first lines in the main BART script

## 0. Load packages ####
source(file.path("src", "required.R"))

## Conditional evaluations
# subset_qc <- TRUE # subset to QC data (optional)

## 1. Load data ####

# Select raster files
env_files <- list(here("data", "raster", "spa_stack.tif"), here("data", "raster", "env_stack.tif"))
spe_files <- list(here("data", "raster", "spa_stack.tif"), here("data", "raster", "distributions_raw.tif"))

# Select QC data (if subset_qc option correctly set)
if (exists("subset_qc") && isTRUE(subset_qc)) {
  message("Subsetting to QC data")
  env_files <- list(here("data", "raster", "spa_stack_qc.tif"), here("data", "raster", "env_stack_qc.tif"))
  spe_files <- list(here("data", "raster", "spa_stack_qc.tif"), here("data", "raster", "distributions_raw_qc.tif"))
}

# Load rasters as stack
(env_stack <- stack(env_files))
(spe_stack <- stack(spe_files))

# Rename variables
names(env_stack) <- c("site", "lon", "lat", paste0("wc", 1:19), paste0("lc", 1:10))
names(spe_stack) <- c("site", "lon", "lat", paste0("sp", 1:(nlayers(spe_stack)-3)))

# Convert to tibble
(env_full <- as_tibble(as.data.frame(env_stack)))
(spe_full <- as_tibble(as.data.frame(spe_stack)))

# Reorder rows by site id (same order as SimpleSDMLayers)
(env_full <- arrange(env_full, site))
(spe_full <- arrange(spe_full, site))

# Select sites with observations only & replace NA by zero
spe <- spe_full %>%
  filter(if_any(contains("sp"), ~ !is.na(.x))) %>%
  mutate(across(contains("sp"), ~ replace(., is.na(.), 0)))
env <- filter(env_full, site %in% spe$site)
spe
env

# Remove site with NAs for landcover variables
(inds_withNAs <- unique(unlist(map(env, ~ which(is.na(.x))))))
if (length(inds_withNAs) > 0) {
  message("Removing sites with observations but NA for land cover values")
  spe <- spe[-inds_withNAs,]
  env <- env[-inds_withNAs,]
}

# Remove species without observations
(spe_withoutobs <- names(which(colSums(spe) == 0)))
if (length(spe_withoutobs) > 0) {
  message("Removing ", length(spe_withoutobs), " species without observations")
  spe <- dplyr::select(spe, -all_of(spe_withoutobs))
  spe_full <- dplyr::select(spe_full, -all_of(spe_withoutobs))
}

# Select fewer variables
xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:3,5,7:10)))
vars_stack <- subset(env_stack, xnames)
