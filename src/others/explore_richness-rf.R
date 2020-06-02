#### Species richness RF model ####

## 0. Load packages ####
source(file.path("src", "required.R"))

# Conditional evaluations
subset_qc <- TRUE # subset to QC data (optional)
# save_models <- TRUE # save & overwrite models


## 1. Prepare data ####

# Load data
source(here("src", "02a_training_data-preparation.R"))

# Get richness values
richness_values <- rowSums(spe)

# Separate into training/testing datasets
set.seed(42)
inds_train <- sample(nrow(vars), 0.7*nrow(vars), replace = FALSE)

richness_train <- richness_values[inds_train]
vars_train <- vars[inds_train,]

richness_test <- richness_values[-inds_train]
vars_test <- vars[-inds_train,]


## 2. Train model ####
# Regression model
system.time(
    regress_model <- ranger(
        richness_train ~ ., 
        data = vars_train, 
        importance = "impurity", 
        seed = 42
    )
)
# Classification model
system.time(
    classif_model <- ranger(
        richness_train ~ ., 
        data = vars_train, 
        classification = TRUE, 
        importance = "impurity", 
        seed = 42
    )
)

# Check accuracy
regress_pred <- predict(regress_model, vars_test)$predictions
sum(round(regress_pred) == richness_test)/length(richness_test)
classif_pred <- predict(classif_model, vars_test)$predictions
sum(classif_pred == richness_test)/length(richness_test)


## 3. Full-scale predictions ####

# Remove sites with NA values
inds_na <- map(env_full, ~ which(is.na(.x)))
(inds_na <- sort(unique(unlist(inds_na))))
vars_nona <- vars_full[-inds_na,]

# Make predictions
predictions_full_scale <- predict(classif_model, vars_nona)$predictions

# Add sites with NAs
predictions_full <- matrix(NA, nrow = nrow(vars_full), ncol = 1)
colnames(predictions_full) <- colnames(predictions_full_scale)
predictions_full[-inds_na,] <- predictions_full_scale