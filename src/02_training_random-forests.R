## 0. Load packages ####
library(conflicted)
library(tidyverse)
library(ranger)
library(caret)
library(pbapply)

# Resolve conflicts
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "dplyr")

# Conditional evaluations
subset_qc <- TRUE # subset to QC data (optional)
# create_models <- TRUE # train models
# save_models <- TRUE # save & overwrite models

## 1. Load data ####

message("Loading & preparing data")

# Load data
spa_full <- read_tsv("data/proc/distributions_spa_full.csv")
env_full <- read_tsv("data/proc/distributions_env_full.csv")
spe      <- read_tsv("data/proc/distributions_spe_full.csv") 

# Load QC data (optional)
spa_qc <- read_tsv("data/proc/distributions_spa_qc.csv")

# Prepare data
# subset_qc <- TRUE # subset to QC data (optional)
source("src/02_training_data-preparation.R")

# Separate into training/testing datasets
set.seed(42)
inds_train <- sample(nrow(spe), 0.7*nrow(spe), replace = FALSE)

spe_train <- spe[inds_train,]
spa_train <- spa[inds_train,]
env_train <- env[inds_train,]
vars_train <- vars[inds_train,]

spe_test <- spe[-inds_train,]
spa_test <- spa[-inds_train,]
env_test <- env[-inds_train,]
vars_test <- vars[-inds_train,]

# Remove species without observations in subsets
(inds_withoutobs <- c(which(sapply(spe_train, sum) == 0), which(sapply(spe_test, sum) == 0)))
if (length(inds_withoutobs > 0)) {
    spe_train <- spe_train[, -inds_withoutobs]
    spe_test <- spe_test[, -inds_withoutobs]
}

# Create single species subset
sp <- "sp17"
sp_train <- spe_train[sp]
sp_test <- as.factor(spe_test[[sp]])

#### Ranger ####

# Single species time test
system.time(
    ranger_model <- ranger(
        sp_train ~ ., 
        data = vars_train,
        classification = TRUE, 
        importance = "impurity", 
        seed = 42
    )
)

# Wrap as function
ranger_train <- function(sp, vars, ...) {
    sp_train <- as.factor(sp)
    rf <- ranger(
        sp_train ~ .,
        data = vars,
        classification = TRUE,
        importance = "impurity",
        seed = 42,
        ...
    )
    return(rf)
}

# Multi-species calls
system.time(
    ranger_models <- pblapply(spe_train, function(x) ranger_train(x, vars_train))
)

# View results
ranger_models

# Combine results in dataframe
rf_res <- tibble(
    species = colnames(spe_train),
    OOB = map_dbl(ranger_models, "prediction.error"),
    error_rate_0 = map_dbl(ranger_models, ~ .x$confusion.matrix[3]/sum(.x$confusion.matrix[c(1,3)])),
    error_rate_1 = map_dbl(ranger_models, ~ .x$confusion.matrix[2]/sum(.x$confusion.matrix[c(2,4)]))
)
rf_res

# Plot results
barplot(rf_res$OOB, names.arg = rf_res$species)
hist(rf_res$OOB, breaks=20)
boxplot(rf_res$OOB)

# Plot presence vs absence error rate
barplot(rf_res$error_rate_0, names.arg = rf_res$species)
barplot(rf_res$error_rate_1, names.arg = rf_res$species)

## Export model
save(ranger_models, file = "data/proc/rf_models.RData")

## Test ranger predictions
# Make predictions
system.time(ranger_tests <- map(ranger_models, ~ predict(.x, vars_test)))
ranger_tests

# Extract confusion matrix and results
ranger_confusion <- map2(
    ranger_tests, spe_test,
    ~ confusionMatrix(
        data = .x$predictions,
        reference = as.factor(.y),
        positive = "1"
    )
)
rf_test_res <- tibble(
    species = colnames(spe_test),
    test_OOB = map_dbl(ranger_confusion, ~ 1 - .x$overall["Accuracy"]),
    test_error_rate_0 = map_dbl(ranger_confusion, ~ 1 - .x$byClass["Specificity"]),
    test_error_rate_1 = map_dbl(ranger_confusion, ~ 1 - .x$byClass["Sensitivity"])
)
rf_test_res

# Compare training & validation errors
(comparison <- (select(rf_res, -species) - select(rf_test_res, -species)) * 100) %>% as_tibble()
summary(comparison)
hist(comparison$OOB) # mostly < 1
hist(comparison$error_rate_0, breaks = 10) # mostly < 1
hist(comparison$error_rate_1, breaks = 20) # some extremes
cor(rf_res[,-1], rf_test_res[,-1]) # correlations of 0.999, 0.999 and 0.982
