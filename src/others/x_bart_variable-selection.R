#### Variable selection ####
# Also see x_bart_single-species.R for a complete single-species attempt

## 0. Load data ####

# Conditional evaluations
subset_qc <- TRUE # subset to QC data (optional)

# Load data
source(file.path("src", "others", "x_training_data.R"))

# Check prepared data
env_full # environmental data for all sites
spe_full # occurrence data for all sites
env # environmental data for sites with observations
spe # occurrence data for sites with observations
env_stack # environmental data layers
spe_stack # species distribution layers
xnames # variables selected for analyses
vars_stack # layers of selected variables

## 1. Variable Selection ###

warning(
  "Testing variable selection for 3 species only.
   This takes 12 min per species for Quebec and 2.5 hours at full scale."
)

# Select fewer species
sort(colSums(spe)/nrow(spe), decreasing = TRUE)
spe_sel <- c("sp1", "sp6", "sp17")
vars_sel <- map(spe_sel, ~ NULL)
names(vars_sel) <- spe_sel

# Run variable selection
tic("total")
for(sp in spe_sel){
  tic(sp)
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
  toc()
}
toc()
vars_sel
