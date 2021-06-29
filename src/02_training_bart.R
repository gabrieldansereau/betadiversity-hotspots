## 0. Load packages ####
source(file.path("src", "required.R"))

## Conditional evaluations
# subset_qc <- TRUE # subset to QC data (optional)
# create_models <- TRUE # train models
# save_models <- TRUE # save & overwrite models
# save_predictions <- TRUE # save & overwrite predictions

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
plot(vars_stack)

## 2. Prepare model functions & outputs ####

message("Preparing model functions & outputs")

# Split species in groups (every group will run on a different core)
# max_procs <- availableCores() - 1
max_procs <- 7 # max number with fits in memory in my case
spe_num <- which(startsWith(names(spe), "sp"))
spe_splits <- split(spe_num, ceiling((spe_num - min(spe_num) + 1)/max_procs))
spe_groups <- map(spe_splits, ~ select(spe, all_of(.)))

# Function to run BART in parallel
bart_parallel <- function(x.train, y.train, ...) {
  # BART SDM
  sdm <- bart(
    y.train = y.train,
    x.train = x.train,
    keeptrees = TRUE,
    verbose = FALSE,
    ...
  )
  # Touch state so that saving will work
  invisible(sdm$fit$state)
  return(sdm)
}

# Prepare global sets
# Models (not currently used)
sdms_list <- list()
predictions_list <- list()
# Summary results for every species
results_global <- tibble(
  spe = names(select(spe, contains("sp"))),
  auc = NA_real_,
  threshold = NA_real_,
  tss = NA_real_,
  type_I = NA_real_,
  type_II = NA_real_
)
# Variable importance for every variable & species
varimps_global <- spe_full %>%
  slice(1:length(xnames)) %>%
  mutate(across(everything(), ~ replace(.x, !is.na(.), NA_real_))) %>%
  mutate(vars = xnames) %>%
  select(vars, everything(), -c(site, lon, lat))
# Predictions assembled as tibbles
pred_df_global <- spe_full %>%
  select(-c(site, lon, lat)) %>%
  mutate(across(everything(), ~ replace(.x, !is.na(.), NA_real_)))
lower_df_global <- pred_df_global
upper_df_global <- pred_df_global
pres_df_global <- pred_df_global

## 3. Run models ####

# The BART model objects take a lot of memory at full scale.
# Model training & predictions are therefore done in the same (huge) loop.
# The loop runs for groups of species which fit in memory at the same time.
# The training & predicting steps are run in parallel.
# Models are exported to .RData files in case of reuse.
# Predictions are exported to the predictions tibbles (and later exported to CSV).

message("Training multi-species models & predicting distributions")

# Run loop for each group
system.time(
for (gp in seq_along(spe_groups)) {

  message(paste0("Training multi-species group (", gp, "/", length(spe_groups)), ")")

  ## 3.1 Create BART models ####

  # Set file paths for .RData
  modelname <- ifelse(exists("subset_qc") && isTRUE(subset_qc), paste0("bart_models_qc", gp, ".RData"), paste0("bart_models", gp, ".RData"))
  filepath <- here("data", "rdata", modelname)
  # Create models (and optionally save to .RData) or load from existing RData
  # create_models <- TRUE
  if (exists("create_models") && isTRUE(create_models)){
    # Run models in parallel
    message("Creating models in parallel: ", modelname)
    set.seed(42)
    system.time(
      sdms <- future_map(
        spe_groups[[gp]],
        ~ bart_parallel(
          y.train = .x,
          x.train = env[, xnames]
        ),
        .progress = TRUE,
        .options = furrr_options(seed = TRUE) # disables warning about seed
      )
    ) # ~ 4 min in parallel

    # Export results
    # save_models <- TRUE
    if (exists("save_models") && isTRUE(save_models)) {
      message("Saving models to RData: ", modelname)
      save(sdms, file = filepath)
    }
  } else {
    # Load models from files
    message("Loading models from RData: ", modelname)
    load(filepath)
  }

  ## 3.2 Extract summary data ####

  # Extract summary statistics
  summaries <-  future_map(sdms, summary_inner)

  # Organize as tibble
  results <- summaries %>%
    map_df(`[`, c("auc", "threshold", "tss", "type_I", "type_II")) %>%
    mutate(spe = names(summaries)) %>%
    select(spe, everything())
  print(results, n = Inf)
  summary(results)

  # Extract variable importance
  varimps <- map(sdms, varimp) %>%
    map_df("varimps") %>%
    mutate(vars = xnames) %>%
    select(vars, everything())
  varimps
  varimps %>%
    pivot_longer(-vars, names_to = "spe", values_to = "varimp") %>%
    pivot_wider(spe, names_from = "vars", values_from = "varimp")
  varimps %>%
    transmute(
      vars = vars,
      mean = rowMeans(select(., -vars))
    ) %>%
    arrange(desc(mean))

  ## 3.3 Predict species distributions ####

  message("Predicting species distributions: ", modelname)

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
  ) # ~ 3 min in parallel

  # Collect predictions
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

  # Convert to presence-absence based on recommended threshold per species
  pres_df <- map2_df(
    pred_df, results$threshold,
    function(pred, thresh) ifelse(pred > thresh, 1, 0)
  )
  pres_df

  ## 3.4 Export results ####

  # Models & predictions objects (not used for now)
  # sdms_list[[gp]] <- sdms
  # predictions_list[[gp]] <- predictions
  # Summary statistics
  spe_names <- names(spe_groups[[gp]])
  results_global[results_global$spe %in% spe_names,] <- results
  varimps_global[names(varimps_global) %in% spe_names] <- select(varimps, -vars)
  # Prediction tibbles
  pred_df_global[names(pred_df_global) %in% spe_names] <- pred_df
  lower_df_global[names(lower_df_global) %in% spe_names] <- lower_df
  upper_df_global[names(upper_df_global) %in% spe_names] <- upper_df
  pres_df_global[names(pres_df_global) %in% spe_names] <- pres_df
}
)

## Check results
# Summary statistics
results_global
# Variable importance
varimps_global
# Predictions (posterior mean)
pred_df_global
# Lower CI predictions
lower_df_global
# Upper CI predictions
upper_df_global
# Presence-absence predictions
pres_df_global

# Export to CSV
# save_predictions <- TRUE
if (exists("save_predictions") && isTRUE(save_predictions)) {
  write_tsv(results_global,  here("data", "proc", "bart_summaries.csv"))
  write_tsv(varimps_global,  here("data", "proc", "bart_varimps.csv"))
  write_tsv(pred_df_global,  here("data", "proc", "bart_predictions_prob.csv"))
  write_tsv(lower_df_global, here("data", "proc", "bart_predictions_lower.csv"))
  write_tsv(upper_df_global, here("data", "proc", "bart_predictions_upper.csv"))
  write_tsv(pres_df_global,  here("data", "proc", "bart_predictions_pres.csv"))
}

## 4. Visualize results ####

# Empty canvas
pred_plot <- ggplot(env_full, aes(lon, lat)) +
  scale_fill_viridis_c(na.value = "white", "Value") +
  coord_quickmap() +
  theme_minimal()

# Plot predictions
sp_no <- 2
pred_plot + geom_raster(aes(fill = pred_df_global[[sp_no]])) + ggtitle("Posterior mean")
pred_plot + geom_raster(aes(fill = lower_df_global[[sp_no]])) + ggtitle("Lower CI")
pred_plot + geom_raster(aes(fill = upper_df_global[[sp_no]])) + ggtitle("Upper CI")
pred_plot + geom_raster(aes(fill = upper_df_global[[sp_no]] - lower_df_global[[sp_no]])) + ggtitle("CI difference")
pred_plot + geom_raster(aes(fill = pres_df_global[[sp_no]])) + ggtitle("Presence-absence")
