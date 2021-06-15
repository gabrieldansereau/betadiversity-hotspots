## 0. Load packages ####
source(file.path("src", "required.R"))

# Conditional evaluations
# subset_qc <- TRUE # subset to QC data (optional)
# create_models <- TRUE # train models
# save_models <- TRUE # save & overwrite models


## 1. Load data ####

# Raster data
(spa_stack <- stack("./data/proc/spa_stack.tif"))
(env_stack <- stack("./data/proc/env_stack.tif"))

names(spa_stack) <- c("site", "lon", "lat")
names(env_stack) <- c(paste0("wc", 1:19), paste0("lc", 1:10))

(spa_full <- as_tibble(as.data.frame(spa_stack)))
(env_full <- as_tibble(as.data.frame(env_stack)))

(vars_full2 <- bind_cols(spa_full, env_full))

vars_full2 <- vars_full2 %>% 
    arrange(site)

# Load data
source(here("src", "02a_training_data-preparation.R"))

vars_full
vars_full2
vars_full == vars_full2
all(vars_full == vars_full2, na.rm = TRUE)
all(select(vars_full, site:lat) == select(vars_full2, site:lat))
drop_na(vars_full) == drop_na(vars_full2)
(idcols <- map_lgl(1:ncol(vars_full), ~ all(vars_full[,.x] == vars_full2[,.x], na.rm = TRUE)))
vars_full[,-which(idcols)]
test <- bind_cols(select(vars_full, wc1), select(vars_full2, wc1)) %>% 
    drop_na()
names(test) <- c("old", "new")
options(pillar.sigfig = 10)
test %>% 
    mutate(eq = old - new) %>% 
    summary()

all.equal(test$old, test$new, tol=1e-7)
all.equal(vars_full, vars_full2, tol=1e-7, na.rm = TRUE)
all.equal(vars_full, vars_full2, tol=1e-7, na.rm = TRUE)
(idcols <- sapply(1:ncol(vars_full), function(x) all.equal(vars_full[,x], vars_full2[,x], tol=1e-7)))
all(sapply(1:ncol(vars_full), function(x) all.equal(vars_full[,x], vars_full2[,x], tol=1e-7)))
# All approximately equal

# Select fewer variables
xnames <- c(paste0("wc", c(1, 2, 5, 6, 12, 13, 14, 15)), paste0("lc", c(1:3,5,7:10)))

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
plot(vars_stack)


## 3. Multi-species model training ####

message("Training multi-species models")

# Split species in groups
# max_procs <- availableCores() - 1
max_procs <- 7
spe_num <- 1:ncol(spe)
spe_splits <- split(spe_num, ceiling(spe_num/max_procs))
spe_groups <- map(spe_splits, ~ select(spe, all_of(.)))

# Function to run BART in parallel
bart_parallel <- function(x.train, y.train, ...) {
    # BART SDM
    sdm <- bart(
        y.train = y.train,
        x.train = x.train,
        keeptrees = TRUE,
        verbose = FALSE,
        ...
    )
    # Touch state so that saving will work
    invisible(sdm$fit$state)
    return(sdm)
}

# Prepare global sets
sdms_list <- list()
predictions_list <- list()
results_global <- tibble(
    spe = names(spe),
    auc = NaN,
    threshold = NaN,
    tss = NaN,
    type_I = NaN,
    type_II = NaN
)
varimps_global <- spe_full %>%
    slice(seq_along(xnames)) %>% 
    mutate_all(~ replace(.x, !is.nan(.), NaN)) %>% 
    mutate(vars = xnames) %>% 
    select(vars, everything(), -site)
pred_df_global <- spe_full %>% 
    select(-site) %>% 
    mutate_all(~ replace(.x, !is.nan(.), NaN))
lower_df_global <- pred_df_global
upper_df_global <- pred_df_global
pres_df_global <- pred_df_global


# Run for each group
system.time(
for (gp in seq_along(spe_groups)) {

    message(paste0("Training multi-species group (", gp, "/", length(spe_groups)), ")")

    ## Create BART models
    # create_models <- TRUE
    modelname <- ifelse(exists("subset_qc") && isTrue(subset_qc), paste0("bart_models_qc", gp, ".RData"), paste0("bart_models", gp, ".RData"))
    filepath <- here("data", "rdata", modelname)
    if (exists("create_models") && isTRUE(create_models)){
        # Run models
        message("Creating models in parallel: ", modelname)
        set.seed(42)
        system.time(
            sdms <- future_map(
                spe_groups[[gp]],
                ~ bart_parallel(
                    y.train = .x,
                    x.train = vars[,xnames]
                ),
                .progress = TRUE
            )
        ) # ~ 4 min in parallel
        
        # Export results
        # save_models <- TRUE
        if (exists("save_models") && isTRUE(save_models)) {
            message("Saving models to RData: ", modelname)
            save(sdms, file = filepath)
        }
    } else {
        # Load models from files
        message("Loading models from RData: ", modelname)
        load(filepath)
    }

    # Extract summary statistics
    summaries <-  future_map(sdms, summary_inner)

    # Organize as tibble
    results <- summaries %>% 
        map_df(`[`, c("auc", "threshold", "tss", "type_I", "type_II")) %>% 
        mutate(spe = names(summaries)) %>% 
        select(spe, everything())
    print(results, n = Inf)
    summary(results)

    # Extract variable importance
    varimps <- map(sdms, varimp) %>% 
        map_df("varimps") %>% 
        mutate(vars = xnames) %>% 
        select(vars, everything())
    varimps
    varimps %>% 
        pivot_longer(-vars, names_to = "spe", values_to = "varimp") %>% 
        pivot_wider(spe, names_from = "vars", values_from = "varimp")
    varimps %>% 
        transmute(
            vars = vars,
            mean = rowMeans(select(., -vars))
        ) %>% 
        arrange(desc(mean))


    ## 4. Multi-species predictions ####

    message("Predicting species distributions: ", modelname)

    # Quantile Predictions
    system.time(
        predictions <- future_map(
            sdms,
            # sdms_backup,
            function(x) predict2.bart(
                object = x, 
                x.layers = vars_stack,
                quantiles = c(0.025, 0.975),
                splitby = 20,
                quiet = TRUE
            ),
            .progress = TRUE
        )
    ) # ~ 3 min in parallel

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

    # Presence-absence dataframe
    pres_df <- map2_df(
        pred_df, results$threshold, 
        function(pred, thresh) ifelse(pred > thresh, 1, 0) 
    )
    pres_df

    # Export group results
    # sdms_list[[gp]] <- sdms
    # predictions_list[[gp]] <- predictions
    results_global[spe_splits[[gp]],] <- results
    varimps_global[spe_splits[[gp]]+1] <- varimps[,-1]
    pred_df_global[spe_splits[[gp]]] <- pred_df
    lower_df_global[spe_splits[[gp]]] <- lower_df
    upper_df_global[spe_splits[[gp]]] <- upper_df
    pres_df_global[spe_splits[[gp]]] <- pres_df
}
)

# Check results
results_global
varimps_global
pred_df_global
lower_df_global
upper_df_global
pres_df_global

# Export to CSV
# save_predictions <- TRUE
if (exists("save_predictions") && isTRUE(save_predictions)) {
    write_tsv(results_global,  here("data", "proc", "bart_summaries.csv"))
    write_tsv(varimps_global,  here("data", "proc", "bart_varimps.csv"))
    write_tsv(pred_df_global,  here("data", "proc", "bart_predictions_prob.csv"))
    write_tsv(lower_df_global, here("data", "proc", "bart_predictions_lower.csv"))
    write_tsv(upper_df_global, here("data", "proc", "bart_predictions_upper.csv"))
    write_tsv(pres_df_global,  here("data", "proc", "bart_predictions_pres.csv"))
}

## 5. Visualize results ####

# Empty canvas
pred_plot <- ggplot(spa_full, aes(lon, lat)) +
    scale_fill_viridis_c(na.value = "white", "Value") +
    coord_quickmap() +
    theme_minimal()

# Plot predictions
sp_no <- 2
pred_plot + geom_raster(aes(fill = pred_df_global[[sp_no]])) + ggtitle("Posterior mean")
pred_plot + geom_raster(aes(fill = lower_df_global[[sp_no]])) + ggtitle("Lower CI")
pred_plot + geom_raster(aes(fill = upper_df_global[[sp_no]])) + ggtitle("Upper CI")
pred_plot + geom_raster(aes(fill = upper_df_global[[sp_no]] - lower_df_global[[sp_no]])) + ggtitle("CI difference")
pred_plot + geom_raster(aes(fill = pres_df_global[[sp_no]])) + ggtitle("Presence-absence")


## Variable selection

# Select species
sort(colSums(spe)/nrow(spe), decreasing = TRUE)
spe_sel <- c("sp1", "sp6", "sp17")
vars_sel <- map(spe_sel, ~ NULL)
names(vars_sel) <- spe_sel
# Run variable selection
# variable_selection <- TRUE
if (exists("variable_selection") && isTRUE(variable_selection)) {
    tictoc::tic("total")
    for(sp in spe_sel){
        tictoc::tic(sp)
        set.seed(42)
        message(paste0("Variable selection for ", sp, " (", which(sp == spe_sel), "/", length(spe_sel)), ")")
        # Save plot to png
        png(here("fig", "bart", paste0("x_bart_vars-select_", sp, ".png")))
        step_vars <- variable.step(
            y.data = spe[[sp]], 
            x.data = env[xnames],
            iter = 50
        )
        dev.off()
        # Save variables to list
        vars_sel[[sp]] <- step_vars
        tictoc::toc()
    }
    tictoc::toc()
    vars_sel
}

## 3 species QC scale
# $sp1
#  [1] "wc1"  "wc2"  "wc5"  "wc6"  "wc12" "wc13" "wc15" "lc1"  "lc2"  "lc3"  "lc5"  "lc8" 
# $sp17
#  [1] "wc1"  "wc2"  "wc5"  "wc6"  "wc14" "wc15" "lc2"  "lc3"  "lc5"  "lc8" 
# $sp9
# [1] "wc1"  "wc5"  "wc6"  "wc15" "lc3"  "lc8"  "lc10"

## 3 species full scale
# $sp1
#  [1] "wc1"  "wc2"  "wc5"  "wc6"  "wc13" "wc14" "wc15" "lc2"  "lc3"  "lc5"  "lc7"  "lc8"  "lc9" 
# $sp6
#  [1] "wc1"  "wc2"  "wc5"  "wc6"  "wc12" "wc14" "wc15" "lc2"  "lc3"  "lc5"  "lc7"  "lc8"  "lc9" 
# $sp17
# [1] "wc1"  "wc5"  "wc6"  "wc12" "wc14" "wc15" "lc7"  "lc8"  "lc9"