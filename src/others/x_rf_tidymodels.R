## 0. Load packages ####
source(file.path("src", "required.R"))

# Conditional evaluations
subset_qc <- TRUE # subset to QC data (optional)
# create_models <- TRUE # train models
# save_models <- TRUE # save & overwrite models

## 1. Load data ####

# Load data
source(here("src", "02a_training_data-preparation.R"))

## 2. Tidymodels attempt
# Inspired by https://juliasilge.com/blog/sf-trees-random-tuning/

library(tidymodels)

# Combine occ & env data for 1 species 
sp_df <- spe %>% 
    select(sp = sp17) %>% 
    mutate_at("sp", factor) %>% 
    bind_cols(vars)
sp_df

# Plot distribution
sp_df %>% 
    ggplot(aes(lon, lat, color = sp)) +
    geom_point(size = 0.5, alpha = 0.4) +
    labs(color = NULL)

# Relationships
sp_df %>% 
    dplyr::count(sp)

## Build model

# Split datasets
set.seed(42)
sp_split <- initial_split(sp_df, strata = sp)
sp_train <- training(sp_split)
sp_test <- testing(sp_split)

# Prepare recipe (but no specific steps in this case)
sp_rec  <- recipe(sp ~ ., data = sp_train)
sp_prep <- prep(sp_rec)
juiced <- juice(sp_prep)

# Tuning parameters
tune_spec <- rand_forest(
    mtry = tune(),
    trees = 1000,
    min_n = tune()
) %>% 
    set_mode("classification") %>% 
    set_engine("ranger")

# Prepare workflow
tune_wf <- workflow() %>% 
    add_recipe(sp_rec) %>% 
    add_model(tune_spec)


## Train hyperparameters

# Cross-validation resamples
set.seed(234)
sp_folds <- vfold_cv(sp_train)

# Enable parallel processing
doParallel::registerDoParallel()

# Create grid
set.seed(345)
tune_res <- tune_grid(
    tune_wf,
    resamples = sp_folds,
    grid = 20
)
tune_res

# Look at performance
tune_res %>% 
    collect_metrics() %>% 
    filter(.metric == "roc_auc") %>% 
    select(mean, min_n, mtry) %>% 
    pivot_longer(min_n:mtry,
        values_to = "value",
        names_to = "parameter"
    ) %>% 
    ggplot(aes(value, mean, color = parameter)) +
    geom_point(show.legend = FALSE) +
    facet_wrap(~parameter, scales = "free_x") +
    labs(x = NULL, y = "AUC")

## Refine hyperparameter grid

# # Define regular grid
# rf_grid <- grid_regular(
#   mtry(range = c(10, 30)),
#   min_n(range = c(2, 8)),
#   levels = 5
# )
# rf_grid

# # Tune
# set.seed(456)
# regular_res <- tune_grid(
#   tune_wf,
#   resamples = trees_folds,
#   grid = rf_grid
# )
# regular_res

# Define regular_res if previous code not run
regular_res <- tune_res
# Plot results
regular_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  mutate(min_n = factor(min_n)) %>%
  ggplot(aes(mtry, mean, color = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point() +
  labs(y = "AUC")


## Choosing the best model
best_auc <- select_best(regular_res, "roc_auc")
final_rf <- finalize_model(
    tune_spec,
    best_auc
)

# Check variable importance
library(vip)
final_rf %>% 
    set_engine("ranger", importance = "permutation") %>% 
    fit(
        sp ~ ., 
        data = juice(sp_prep)
    ) %>% 
    vip(geom = "point")

# Build final workflow
final_wf <- workflow() %>% 
    add_recipe(sp_rec) %>% 
    add_model(final_rf)
# Get final results on test set
final_res <- final_wf %>% 
    last_fit(sp_split)
# Get final metrics on test set
final_res %>% 
    collect_metrics()

# Map incorrectly labeled (in test set)
final_res %>% 
    collect_predictions() %>% 
    mutate(
        correct = case_when(
            sp == .pred_class ~ "Correct",
            TRUE ~ "Incorrect"
        )
    ) %>% 
    bind_cols(sp_test %>% select(-sp)) %>% 
    ggplot(aes(lon, lat, color = correct)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(color = NULL) +
    scale_color_manual(values = c("gray80", "darkred"))


## Custom : parallel training for single model (instead of during crossvalidation)

# Define model
cores <- parallel::detectCores()
parallel_rf <- rand_forest(
    mtry = 6,
    trees = 1000,
    min_n = 10
) %>% 
    set_mode("classification") %>% 
    set_engine("ranger", num.threads = cores)

# Attempt to fit
parallel_res <- parallel_rf %>% 
    fit(sp ~ ., data = juice(sp_prep))

# Put in workflow
parallel_wf <- workflow() %>% 
    add_recipe(sp_rec) %>% 
    add_model(parallel_rf)

# Results on test set
parallel_res <- parallel_wf %>% 
    last_fit(sp_split)

parallel_res %>% 
    collect_metrics()

test_predictions <- parallel_res %>% 
    collect_predictions()

test_predictions %>% 
    conf_mat(truth = sp, estimate = .pred_class)

test_predictions %>%
  ggplot() +
  geom_density(aes(x = .pred_1, fill = sp), 
               alpha = 0.5)
# Attempt to predict at full scale
parallel_model <- fit(parallel_wf, sp_train) # option 1a
parallel_model <- fit(parallel_wf, juice(sp_prep)) # option 1b # preferred one if recipe changes data
parallel_model <- fit(parallel_rf, sp ~ ., data = sp_train) # option 2a
parallel_model <- fit(parallel_rf, sp ~ ., data = juice(sp_prep)) # option 2b
predicted <- predict(
    parallel_model, 
    new_data = vars_full %>% na.omit()
)

# Plot predictions
predicted %>% 
    bind_cols(vars_full %>% na.omit()) %>% 
    ggplot(aes(lon, lat, color = .pred_class)) +
    geom_point(size = 0.5, alpha = 0.4) +
    scale_color_viridis_d() +
    labs(color = NULL)