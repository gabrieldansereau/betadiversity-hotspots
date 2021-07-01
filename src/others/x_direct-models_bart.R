#### Direct richness/LCBD BART models ####

## 0. Load data ####

# Conditional evaluations
subset_qc <- TRUE # subset to QC data (optional)

# Load data
source(file.path("src", "others", "x_training_data.R"))

# Check prepared data
env_full # environmental data for all sites
spe_full # occurrence data for all sites
env # environmental data for sites with observations
spe # occurrence data for sites with observations
env_stack # environmental data layers
spe_stack # species distribution layers
xnames # variables selected for analyses
vars_stack # layers of selected variables

## 1. Prepare data ####

# Get richness values
richness_values <- rowSums(select(spe, starts_with("sp")))

# Get LCBD values
# Apply Hellinger transformation
sites_obs <- spe$site
inds_obs  <- which(env_full$site %in% sites_obs)
spe_sp <- select(spe, contains("sp"))
spe_transf <- vegan::decostand(spe_sp, "hel")
# S -> squared deviations from column mean
S <- map_df(spe_transf, ~ (.x - mean(.x))^2)
# SStotal -> total sum of squares
SStotal <- sum(S)
# SSi -> sum of squares for site i
SSi <- rowSums(S)
# LCBD -> local contribution to beta diversity (site i, relative)
lcbd_values <- SSi/SStotal
# Set values relative to maximum
lcbd_values <- lcbd_values/max(lcbd_values)

## 2. Model training ####

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

# Train models
set.seed(42)
tic("Training direct models")
models <- future_map(
  values_df,
  ~ bart_parallel(
    y.train = .x,
    x.train = env[,xnames],
    keeptrees = TRUE
  ),
  .options = furrr_options(seed = TRUE) # disables warning about seed
); toc() # 90 sec.

# Check variable importance
varimp(models[[1]], plot = TRUE)
varimps <-  map_dfr(models, varimp, .id = "value") %>%
  rename(vars = names)

# Compare variable importance in richness vs LCBD model
varimps %>%
  pivot_wider(names_from = "value", values_from = "varimps") %>%
  mutate(diff = richness - lcbd)


ggplot(varimps, aes(vars)) +
  geom_bar(aes(weight = varimps, fill = value)) #, position = "dodge2")

ggplot(varimps, aes(vars, varimps)) +
  geom_point(aes(colour = value))

## 3. Predictions ####

# Predict values
tic("Predicting from direct models")
predictions <- future_map(
  models,
  ~ predict2.bart(
    .x,
    vars_stack,
    quantiles = c(0.025, 0.0975),
    splitby = 20
  )
); toc() # 2 min.

# Plot richness predictions
plot(
  predictions$richness[[1]],
  main = 'Richness predictions - Posterior mean',
  col = viridis(255),
  legend.args = list(text='Richness', side=2, line=1.3),
  box = FALSE,
  axes = FALSE
)
# Negative richness??
plot(
  predictions$richness[[1]] < 0,
  main = 'Richness predictions - Negative values',
  col = viridis(255),
  legend.args=list(text='Negative (true/false)', side=2, line=1.3),
  box = FALSE,
  axes = FALSE
)

# Plot LCBD predictions
plot(
  predictions$lcbd[[1]],
  main = 'LCBD predictions - Posterior mean',
  col = viridis(255),
  zlim = c(0, 1),
  legend.args = list(text='Probability', side=2, line=1.3),
  box = FALSE,
  axes = FALSE
)

## Collect predictions

# Mean predictions
pred_df <- predictions %>%
  map(~ .x$layer.1) %>%
  stack() %>%
  as.data.frame(xy = TRUE) %>%
  as_tibble() %>%
  arrange(x, y) %>%
  rename(longitude = x, latitude = y)
pred_df

# Lower quantiles
lower_df <- predictions %>%
  map(~ .x$layer.2) %>%
  stack() %>%
  as.data.frame(xy = TRUE) %>%
  as_tibble() %>%
  arrange(x, y) %>%
  rename(longitude = x, latitude = y)
lower_df

# Upper quantiles
upper_df <- predictions %>%
  map(~ .x$layer.3) %>%
  stack() %>%
  as.data.frame(xy = TRUE) %>%
  as_tibble() %>%
  arrange(x, y) %>%
  rename(longitude = x, latitude = y)
upper_df

# Combine results in list
results <- list(predictions = pred_df, lower = lower_df, upper = upper_df)
results
