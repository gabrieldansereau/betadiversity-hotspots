#### Single-species BART model ####

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

## 1. Single species BART model ####

# Select 1 species with ~ same number of presences & absences
colSums(select(spe, starts_with("sp")))/nrow(select(spe, starts_with("sp"))) # 17 seems good
sp <- select(spe, sp17) # black-throated blue warbler
sp_full <- select(spe_full, sp17) # black-throated blue warbler
table(sp)

# Train model
set.seed(42)
tic("Training model")
sdm <- bart(
  y.train = sp[[1]],
  x.train = env[,xnames],
  keeptrees = TRUE
); toc() # 5 sec. for QC, 90 sec. at full scale

# Model diagnostics
summary(sdm)
varimp(sdm, plot = TRUE)

# Extract inner values
diagnostics <- summary_inner(sdm)

# Predict species distribution
predictions <- predict2.bart(sdm, vars_stack, quantiles=c(0.025, 0.975), splitby = 20)

## 2. Plot singles species predictions ####

# Plot probability predictions
plot(
  env_stack$wc1,
  main = "Original distribution",
  col = "grey",
  legend = FALSE, box = FALSE, axes = FALSE
)
plot(spe_stack$sp1, add = TRUE, col = viridis(2))
plot(
  predictions[[1]],
  main = 'Probability predictions - Posterior mean',
  col = viridis(255),
  zlim = c(0, 1),
  legend.args=list(text='Probability', side=2, line=1.3),
  box = FALSE, axes = FALSE
)
plot(
  predictions[[2]],
  main = 'Probability predictions - Lower 95% CI bound',
  col = viridis(255),
  zlim = c(0, 1),
  legend.args=list(text='Probability', side=2, line=1.3),
  box = FALSE, axes = FALSE
)
plot(
  predictions[[3]],
  main = 'Probability predictions - Upper 95% CI bound',
  col = viridis(255),
  zlim = c(0, 1),
  legend.args=list(text='Probability', side=2, line=1.3),
  box = FALSE, axes = FALSE
)
plot(
  predictions[[3]] - predictions[[2]],
  main = 'Probability predictions - Credible interval width',
  col = viridis(255),
  zlim = c(0, 1),
  legend.args=list(text='Probability', side=2, line=1.3),
  box = FALSE, axes = FALSE
)

# Plot threshold predictions
threshold <- diagnostics$threshold
plot(
  predictions[[1]] > threshold,
  main = 'Threshold predictions - Posterior mean',
  col = viridis(2),
  box = FALSE, axes = FALSE
)
plot(
  predictions[[2]] > threshold,
  main = 'Threshold predictions - Lower 95% CI bound',
  col = viridis(2),
  box = FALSE, axes = FALSE
)
plot(
  predictions[[3]] > threshold,
  main = 'Threshold predictions - Upper 95% CI bound',
  col = viridis(255),
  box = FALSE, axes = FALSE
)
plot(
  (predictions[[3]] > threshold) - (predictions[[2]] > threshold),
  main = 'Threshold predictions - Upper-lower site difference',
  col = viridis(2),
  box = FALSE, axes = FALSE
)

## 3. Variable selection ####

# Variable importance plot
varimp(sdm, plot = TRUE)

# Diagnostic variable selection plot
tic("Diagnostic variable selection")
# varimp.diag(vars[,xnames], sp[[1]], iter = 50)
toc() # ~ 15 min

# Stepwise variable set reduction
tic("Stepwise variable selection")
step.model <- variable.step(
  x.data = env[,xnames],
  y.data = sp[[1]]
); toc() # 12 min for QC, 2:30:00 for full scale
step.model

## 4. Partial dependence ####

# Partial dependence plot
tic("Partial dependence plot")
embarcadero::partial(
  sdm,
  x.vars = c('wc1'),
  trace = FALSE,
  ci = FALSE,
  smooth = 5,
  equal = TRUE
); toc() # ~ 5 min.

# Spartial dependence raster (spatial partial dependence)
tic("Spartial dependence")
p <- spartial(sdm, vars_stack, x.vars='wc1', equal=TRUE, smooth=5); toc() # ~ 3 min.

# Spartial dependence plot
plot(p)

## 5. Others ####
# Not working, can't find function plot.mcmc

# # Show first iterations
# tic("First iterations")
# plot.mcmc(sdm, vars_stack, iter=5, quiet = TRUE); toc()

# # Show burn-in with fewer trees
# set.seed(42)
# sdm.tiny <- bart(
#   y.train = sp[[1]],
#   x.train = env[,xnames],
#   keeptrees = TRUE,
#   ntree = 5, # 5 tree models
#   nskip = 0
# ) # No burnin
# summary(sdm.tiny)
# plot.mcmc(sdm.tiny, vars_stack, iter=100)

# # Timelapse of tree learning
# library(animation)
# saveGIF(
#   plot.mcmc(sdm, vars_stack, iter=50),
#   movie.name = "Timelapse.gif",
#   # interval = 0.15,
#   ani.width = 800,
#   ani.height = 400
# )
