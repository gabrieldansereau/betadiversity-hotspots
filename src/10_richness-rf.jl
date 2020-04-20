import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load distribution data for all species
@load "data/jld2/raw-distributions.jld2" distributions spenames speindex

## Load matrix Y
@load "data/jld2/raw-Y-matrices.jld2" Y Yobs Ytransf inds_obs inds_notobs

## Richness
richness_raw = calculate_richness(Y, inds_notobs, distributions)

# Visualize
plotSDM(richness_raw, c = :viridis)

## Extract values for model
richness_values = Int64.(richness_raw.grid[inds_obs])

## Train Random Forest
@rput richness_values inds_obs
begin
  R"""
  library(ranger)
  spa <-  read.csv("data/proc/distributions_spa.csv", header=TRUE, sep="\t")
  env <-  read.csv("data/proc/distributions_env.csv", header=TRUE, sep="\t")

  # Remove site with NAs for landcover variables
  (inds_withNAs <- unique(unlist(sapply(env, function(x) which(is.na(x))))))
  if (length(inds_withNAs) > 0) {
    richness_values <- richness_values[-inds_withNAs]
    spa <- spa[-inds_withNAs,]
    env <- env[-inds_withNAs,]
  }

  # Combine environmental variables
  vars <- cbind(env, spa)

  # Separate into training/testing datasets
  set.seed(42)
  inds_train <- sample(nrow(vars), 0.7*nrow(vars), replace = FALSE)

  richness_train <- richness_values[inds_train]
  vars_train <- vars[inds_train,]

  richness_test <- richness_values[-inds_train]
  vars_test <- vars[-inds_train,]

  # Train model
  system.time(regress_model <- ranger(richness_train ~ ., data = vars_train, importance = "impurity", seed = 42))
  system.time(classif_model <- ranger(as.factor(richness_train) ~ ., data = vars_train, importance = "impurity", seed = 42))

  regress_pred <- predict(regress_model, vars_test)$predictions
  sum(round(regress_pred) == richness_test)/length(richness_test)

  classif_pred <- predict(classif_model, vars_test)$predictions
  sum(classif_pred == richness_test)/length(richness_test)

  """
end

begin
  R"""
  predictions <- predict(classif_model, vars)$predictions
  predictions <- as.numeric(levels(predictions))[predictions]
  """
end

@rget predictions inds_withNa

richness_rf = similar(richness_raw)
richness_raw.grid[inds_obs] .= predictions
