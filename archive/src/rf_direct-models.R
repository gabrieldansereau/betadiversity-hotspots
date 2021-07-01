#### Species richness RF model ####

## 0. Load packages ####
source(file.path("src", "required.R"))

# Conditional evaluations
# subset_qc <- TRUE # subset to QC data (optional)
# save_models <- TRUE # save & overwrite models


## 1. Prepare data ####

# Load data
source(here("src", "02a_training_data-preparation.R"))

## Get richness values
richness_values <- rowSums(spe)

## Get LCBD values
# Apply Hellinger transformation
spe_transf <- vegan::decostand(spe_full[inds_obs,-1], "hel")
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
# Remove sites with NaN for landcover variables
if (length(inds_withNAs) > 0) {
    lcbd_values <- lcbd_values[-inds_withNAs]
}

# Group richness & lcbd
values_df <- tibble(
    richness = richness_values,
    lcbd = lcbd_values
)

# Separate into training/testing datasets
set.seed(42)
inds_train <- sample(nrow(vars), 0.7*nrow(vars), replace = FALSE)

values_train <- values_df[inds_train,]
vars_train <- vars[inds_train,]

values_test <- values_df[-inds_train,]
vars_test <- vars[-inds_train,]


## 2. Train model ####
# Regression models
system.time(
    models <- map(
        values_train,
        ~ ranger(
            .x ~ ., 
            data = vars_train, 
            importance = "impurity", 
            seed = 42
        )
    )
)

# Check accuracy
pred_test <- map_df(models, ~ predict(.x, vars_test)$predictions)
sum(round(pred_test$richness) == values_test$richness)/length(values_test$richness)


## 3. Full-scale predictions ####

# Remove sites with NA values
inds_na <- map(env_full, ~ which(is.na(.x)))
(inds_na <- sort(unique(unlist(inds_na))))
vars_nona <- vars_full[-inds_na,]

# Make predictions
predictions_nona <- map_df(models, ~ predict(.x, vars_nona)$predictions)

# Add sites with NAs
predictions <- tibble(
    richness = as.numeric(NA),
    lcbd = as.numeric(NA),
    .rows = nrow(vars_full)
)
predictions[-inds_na,] <- predictions_nona

# Put final result in list
results <- list(predictions = predictions)