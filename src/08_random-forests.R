library(party)
library(caret)
library(randomForest)
library(parallel)
library(ranger)
library(pbapply)

#### 0. Load data ####
# Load QC data
spe <- read.csv("data/proc/distributions_spe_qc.csv", header=TRUE, sep="\t")
spa <- read.csv("data/proc/distributions_spa_qc.csv", header=TRUE, sep="\t")
env <- read.csv("data/proc/distributions_env_qc.csv", header=TRUE, sep="\t")

# Remove site with NAs for landcover variables
(inds_withNAs <- unique(unlist(sapply(env, function(x) which(is.na(x))))))
if (length(inds_withNAs) > 0) {
    spe <- spe[-inds_withNAs,]
    spa <- spa[-inds_withNAs,]
    env <- env[-inds_withNAs,]
}

# Combine environmental variables
vars <- cbind(env, spa)

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
sp <- "sp1"
sp_train <- as.factor(spe_train[,sp])
sp_test <- as.factor(spe_test[,sp])

#### 1. Party ####

rf <- cforest(as.factor(spe$sp20) ~ .,
              data = vars,
              control = cforest_unbiased(mtry = 2, ntree = 100))

(vars_imp <- varimp(rf, conditional = T))
(vars_imp <- varimp(rf))

sp1_pred <- predict(rf, OOB=TRUE)

table(spe$sp1, sp1_pred)

#### 2. randomForest ####

set.seed(42)
system.time(
    model1 <- randomForest(sp_train ~ .,
                            data = vars_train,
                            importance = TRUE)
)
model1
# Predict test set
pred_test <- predict(model1, vars_test, type = "class")
# Checking classification accuracy
table(sp_test, pred_test)
mean(sp_test == pred_test)*100

# Check variable importance
importance(model1)
varImpPlot(model1)

 # Select number of trees
set.seed(42)
model2 <- randomForest(sp_train ~ .,
                       data = vars_train,
                       importance = TRUE,
                       ntree = 1000)
set.seed(42)
model3 <- randomForest(sp_train ~ .,
                       data = vars_train,
                       importance = TRUE,
                       ntree = 2000)
par(mfrow=c(2,2))
plot(model1$err.rate[,"OOB"], type = "l", xlab = "ntree", ylab = "OOB error rate")
plot(model2$err.rate[,"OOB"], type = "l", xlab = "ntree", ylab = "OOB error rate")
plot(model3$err.rate[,"OOB"], type = "l", xlab = "ntree", ylab = "OOB error rate")
par(mfrow=c(1,1))

## Find optimal mtry value
set.seed(42)
mtry <- tuneRF(vars_train, sp_train, ntreeTry=500,
               stepFactor=1.5,improve=0.01, trace=TRUE, plot=TRUE)
# Best mtry = 5, same as default

## Wrap as function
rf_train <- function(sp, vars, ...) {
    set.seed(42)
    sp_train <- as.factor(sp)
    rf <- randomForest(sp_train ~ .,
                       data = vars,
                       importance = TRUE,
                       ...)
  return(rf)
}

system.time(
    rf_models <- mclapply(spe_train, function(x) rf_train(x, vars_train), mc.cores = 12)
    )

rf_res <- data.frame(species = colnames(spe_train),
                     OOB = sapply(rf_models, function(x) 1 - sum(x$y == x$predicted)/length(x$y)),
                     error_rate_0 = sapply(rf_models, function(x) median(x$err.rate[,"0"])),
                     error_rate_1 = sapply(rf_models, function(x) median(x$err.rate[,"1"]))
                     )
rf_res

barplot(rf_res$OOB, names.arg = rf_res$species)
hist(rf_res$OOB, breaks=20)
boxplot(rf_res$OOB)

# Test on testing subset
system.time(rf_tests <- lapply(rf_models, function(x) predict(x, vars_test, type = "class")))
rf_tests

rf_res$test_error_rate <- mapply(function(x,y) 1 - confusionMatrix(as.factor(x), y)$overall["Accuracy"], spe_test, rf_tests)
barplot(rf_res$test_error_rate, names.arg = rf_res$species)

#### 3. Ranger ####
# Time models
system.time(rf_model <- randomForest(sp_train ~ ., data = vars_train, importance = TRUE))
system.time(ranger_model <- ranger(sp_train ~ ., data = vars_train, importance = "impurity", seed = 42))

## Wrap as function
ranger_train <- function(sp, vars, ...) {
    sp_train <- as.factor(sp)
    rf <- ranger(sp_train ~ .,
                 data = vars,
                 importance = "impurity",
                 seed = 42,
                 ...)
    return(rf)
}

# Multi-species calls
system.time(
    rf_models <- lapply(spe_train[,1:10], function(x) rf_train(x, vars_train))
)
system.time(
    ranger_models <- pblapply(spe_train, function(x) ranger_train(x, vars_train))
)

# Parallized calls
system.time(
    rf_models <- mclapply(spe_train[,], function(x) rf_train(x, vars_train), mc.cores = 12)
)
system.time(
    ranger_models <- mclapply(spe_train[,], function(x) ranger_train(x, vars_train), mc.cores = 12)
)

# View results
rf_models
ranger_models

# Combine results in dataframe
rf_res <- data.frame(species = colnames(spe_train),
                     OOB = sapply(ranger_models, function(x) x$prediction.error),
                     error_rate_0 = sapply(ranger_models, function(x) x$confusion.matrix[3]/sum(x$confusion.matrix[c(1,3)])),
                     error_rate_1 = sapply(ranger_models, function(x) x$confusion.matrix[2]/sum(x$confusion.matrix[c(2,4)]))
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
system.time(ranger_tests <- lapply(ranger_models, function(x) predict(x, vars_test)))
ranger_tests

# Extract confusion matrix and results
ranger_confusion <- mapply(function(x,y) confusionMatrix(data = x$predictions,
                                           reference = as.factor(y),
                                           positive = "1"),
                           ranger_tests, spe_test, SIMPLIFY = FALSE)
rf_test_res <- data.frame(species = colnames(spe_test),
                          test_OOB = sapply(ranger_confusion, function(x) 1 - x$overall["Accuracy"]),
                          test_error_rate_0 = sapply(ranger_confusion, function(x) 1 - x$byClass["Specificity"]),
                          test_error_rate_1 = sapply(ranger_confusion, function(x) 1 - x$byClass["Sensitivity"]))
head(rf_test_res)

# Compare training & validation errors
(comparison <- (rf_res - rf_test_res) * 100)
summary(comparison)
hist(comparison$OOB) # mostly < 1
hist(comparison$error_rate_0, breaks = 10) # mostly < 1
hist(comparison$error_rate_1, breaks = 20) # some extremes
cor(rf_res[,-1], rf_test_res[,-1]) # correlations of 0.999, 0.999 and 0.982

## Training on full dataset
system.time(
    full_models <- pblapply(spe, function(x) ranger_train(x, vars))
)
full_models

# Combine results in dataframe
rf_full_res <- data.frame(species = colnames(spe),
                     OOB = sapply(full_models, function(x) x$prediction.error),
                     error_rate_0 = sapply(full_models, function(x) x$confusion.matrix[3]/sum(x$confusion.matrix[c(1,3)])),
                     error_rate_1 = sapply(full_models, function(x) x$confusion.matrix[2]/sum(x$confusion.matrix[c(2,4)]))
)

# Compare full training with 70% training
(full_comparison <- (rf_res - rf_full_res[-60,]) * 100)
summary(full_comparison) # not much prediction gain
hist(full_comparison$OOB, breaks = 20) # mostly < 1
hist(full_comparison$error_rate_0, breaks = 20) # mostly < 1
hist(full_comparison$error_rate_1, breaks = 20) # some extremes
cor(rf_res[,-1], rf_full_res[-60,-1])
