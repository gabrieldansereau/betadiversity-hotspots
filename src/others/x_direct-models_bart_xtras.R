#### Direct richness/LCBD BART models - Extras ####

## x0. Prepare data ####

# Conditional evaluations
subset_qc <- TRUE # subset to QC data (optional)

# Run direct model scrip
source(file.path("src", "others", "x_direct-models_bart.R"))

## x1. Variable selection ####

# Stepwise variable reduction
set.seed(42)
tic("Stepwise variable reduction")
step_vars <- variable.step(
  y.data = values_df[[1]],
  x.data = env[xnames],
  iter = 50
); toc() # 7.5 min

# Save variable selection results
step_xvars20_qc <- c("wc1", "wc2", "wc5", "wc6", "wc15", "lc2", "lc3", "lc8")
step_xvars50_full <- c("wc5", "wc6", "wc15")
step_xvars50_full_alternative <- c("wc1", "wc2", "wc5", "wc6", "wc15", "lc3", "lc8", "lc7")

# Train new model
set.seed(42)
sdm_step <- bart(y.train = values_df[[1]], x.train = env[step_vars], keeptrees = TRUE)
summary(sdm_step)
summary(sdm)

# One step selection
set.seed(42)
tic("One step selection")
full_step <- bart.step(
  y.data = values_df[[1]],
  x.data = env[xnames],
  iter.step = 10,
  full = FALSE
); toc() # Error with ROCR on regression model, only supports classification

## Custom function for selection

# Create function to chain variable selection & re-running model
bart.step_custom <- function(x.data, y.data, iter = 50, ...) {
  # Variable selection
  step_vars <- variable.step(y.data = y.data, x.data = x.data, iter = iter)
  # Model on selected variables
  step_model <- bart(y.train = y.data, x.train = x.data[step_vars], keeptrees = TRUE, ...)
  # Touch state so that saving will work
  invisible(step_model$fit$state)
  return(step_model)
}

# Apply on richness
set.seed(42)
tic("Custom function - Richness only")
full_step <- bart.step_custom(
  y.data = values_df[[1]],
  x.data = env[xnames],
  iter = 10
); toc() # 1.5 min for QC

# Apply on richness & LCBD
set.seed(42)
tic("Custom function - Richness & LCBD")
step_models <- future_map(
  values_df,
  ~ bart.step_custom(
    y.data = .x,
    x.data = env[xnames],
    iter = 10
  ),
  .options = furrr_options(seed = TRUE) # disables warning about seed
); toc() # 1.5 min

# Get results
step_vars <- attr(step_models[[1]]$fit$data@x, "term.labels")

## x. Partial dependence plots ####

# Richness
embarcadero::partial(models$richness, x.vars = "wc1", trace = FALSE)
# Not working with continuous data ?? Axis is only between 0 & 1??

# LCBD
embarcadero::partial(models$lcbd, x.vars = "wc1", trace = FALSE)
# Ok, works only for values between 0-1 ?? Is there a way to change y-axis ??

# In parallel
xnames_valid <- names(select(env, all_of(xnames), -lc6))
tic("Partial dependence")
partdep <- future_map(
  xnames_valid,
  ~ embarcadero::partial(
    models$lcbd,
    x.vars = .x,
    equal = TRUE,
    smooth = 1, # less smoothing is faster
    ci = TRUE,
    trace = FALSE # bugs if not set to FALSE in my case
  ),
); toc() # 1.5 min for QC data
names(partdep) <- xnames_valid

# Why is it a list of lists? Let's have a list of plots
partdep <- map(partdep, 1)
partdep$wc1

# Combine 2 plots
partdep$wc1 + partdep$wc12

# Combine all plots
wrap_plots(plotlist = partdep)

# Spartial dependence plots
spartdep <- map(
  xnames_valid,
  ~ spartial(models$lcbd, vars_stack[[xnames]], x.vars=.x, equal=TRUE)
) # better not to run this in parallel or it overflows memory
names(spartdep) <- xnames_valid
spartdep_stack <- stack(spartdep)
plot(spartdep_stack)
plot(spartdep[[1]])
