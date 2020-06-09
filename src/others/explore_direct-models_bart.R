#### Species richness BART model ####

## 0. Load packages ####
source(file.path("src", "required.R"))

# Conditional evaluations
# subset_qc <- TRUE # subset to QC data (optional)
# save_models <- TRUE # save & overwrite models

## 1. Prepare data ####

# Load data
source(here("src", "02a_training_data-preparation.R"))

# Get richness values
richness_values <- rowSums(spe)

# Get LCBD values
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

# Select variables
xnames <- select(vars, -c(lat, lon)) %>% names()
xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:10)))


## 2. Create layers ####

message("Creating layers")

# Create raster layers
vars_layers <- map(
    vars_full[,xnames], 
    ~ df_to_layer(.x, lons = vars_full$lon, lats = vars_full$lat)
)
wc_layer <- vars_layers$wc1

# Stack variables layers
vars_stack <- stack(vars_layers, names = xnames)


## 3. Model training ####

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

# Train model
set.seed(42)
system.time(
    models <- future_map(
        values_df,
        ~ bart_parallel(
            y.train = .x,
            x.train = vars[,xnames],
            keeptrees = TRUE
        )
    )
) # 90 sec.
varimp(models[[1]], plot = TRUE)
varimps <-  map_dfr(models, varimp, .id = "value") %>% 
    rename(vars = names)
varimps %>% 
    pivot_wider(names_from = "value", values_from = "varimps") %>% 
    mutate(diff = richness - lcbd)
ggplot(varimps, aes(vars)) + 
    geom_bar(aes(weight = varimps, fill = value)) #, position = "dodge2")

ggplot(varimps, aes(vars, varimps)) + 
    geom_boxplot(aes(colour = value))


## 4. Predictions ####
system.time(
    predictions <- future_map(
        models,
        ~ predict2.bart(
            .x,
            vars_stack,
            quantiles = c(0.025, 0.0975),
            splitby = 20
        )
    )
) # 2 min.

# Plot richness
plot(
    predictions$richness[[1]], 
    main = 'Probability predictions - Posterior mean', 
    col = viridis(255),
    # zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)
# Negative richness??
plot(
    predictions$richness[[1]] < 1, 
    main = 'Probability predictions - Posterior mean', 
    col = viridis(255),
    # zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)
plot(
    predictions$lcbd[[1]], 
    main = 'Probability predictions - Posterior mean', 
    col = viridis(255),
    zlim = c(0, 1),
    legend.args=list(text='Probability', side=2, line=1.3),
    box = FALSE, axes = FALSE
)

# Predictions
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

# Combine results in list
results <- list(predictions = pred_df, lower = lower_df, upper = upper_df)


## x. Variable selection ####

# Stepwise variable reduction
set.seed(42)
system.time(
    step_vars <- variable.step(
        y.data = values_df[[1]], 
        x.data = env[xnames], 
        iter = 50
    )
) # 3.5 min with 20 trees
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
system.time(
    full_step <- bart.step(
        y.data = values_df[[1]], 
        x.data = env[xnames],
        iter.step = 10,
        full = FALSE
    )
) # Error with ROCR on regression model

# Custom function
bart.step_custom <- function(x.data, y.data, iter = 50, ...) {
    # Variable selection
    step_vars <- variable.step(y.data = y.data, x.data = x.data, iter = iter)
    # Model on selected variables
    step_model <- bart(y.train = y.data, x.train = x.data[step_vars], keeptrees = TRUE, ...)
    # Touch state so that saving will work
    invisible(step_model$fit$state)
    return(step_model)
}
set.seed(42)
system.time(
    full_step <- bart.step_custom(
        y.data = values_df[[1]], 
        x.data = env[xnames], 
        iter = 10
    )
)
set.seed(42)
system.time(
    step_models <- future_map(
        values_df,
        ~ bart.step_custom(
            y.data = .x, 
            x.data = env[xnames], 
            iter = 10
        )
    )
)
step_vars <- attr(step_models[[1]]$fit$data@x, "term.labels")


## x. Partial dependence plots ####

embarcadero::partial(models$richness, x.vars = "wc1", trace = FALSE)
# Not working with continuous data ??
embarcadero::partial(models$lcbd, x.vars = "wc1", trace = FALSE)
# Ok, works only for values between 0-1 ?? Is there a way to change y-axis ??
system.time(
    partdep <- future_map(
        names(select(env, all_of(xnames), -lc6)),
        ~ embarcadero::partial(
            models$lcbd, 
            x.vars = .x,
            equal = TRUE,
            smooth = 1, # less smoothing is faster
            ci = TRUE,
            trace = FALSE # bugs if not set to FALSE in my case
        ),
        .progress = TRUE
    )
) # 2.5 min for QC data, peak at 25 GB Ram...
names(partdep) <- names(select(env, all_of(xnames), -lc6))
partdep[[11]]
par(mfrow = c(2,2))
partdep[[1]]
# How to plot these side by side ??
plot(predictions[[1]]$layer.1)
# Oh now it's working ??
par(mfrow = c(1,1))

spartdep <- spartial(models$lcbd, vars_stack[[xnames]], x.vars="wc1", equal=TRUE)
plot(spartdep)
