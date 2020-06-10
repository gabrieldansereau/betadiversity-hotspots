source(file.path("src", "others", "explore_direct-models_bart.R"))

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
xnames_valid <- names(select(env, all_of(xnames), -lc6))
system.time(
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
        .progress = TRUE
    )
) # 2.5 min for QC data, peak at 25 GB Ram...
names(partdep) <- xnames_valid
# Why is it a list of lists? Let's have a list of plots
partdep <- map(partdep, 1)
partdep$wc1

# Combine 2 plots
gridExtra::grid.arrange(partdep_nolist$wc1, partdep_nolist$wc12, ncol = 2)
cowplot::plot_grid(partdep$wc1, partdep$wc12)
partdep$wc1 + partdep$wc12
# Combine all plots
gridExtra::grid.arrange(plotlist = partdep) # not working
cowplot::plot_grid(plotlist = partdep)
patchwork::wrap_plots(plotlist = partdep)

# Spartial dependence plots
spartdep <- future_map(
    xnames_valid,
    ~ spartial(models$lcbd, vars_stack[[xnames]], x.vars=.x, equal=TRUE),
    .progress = TRUE
)
names(spartdep) <- xnames_valid
spartdep_stack <- stack(spartdep)
plot(spartdep_stack)
plot(spartdep[[1]])
