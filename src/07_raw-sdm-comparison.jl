import Pkg; Pkg.activate(".")
include("required.jl")

# Run analysis on raw data
outcome = "raw"
include("04_analysis.jl")

raw = (Y = Y,
       richness = richness,
       richness_plot = richness_plot,
       lcbd = lcbd,
       beta_total = beta_total,
       lcbd_plot = lcbdtr_plot,
       rel_plot = rel2d_plot
       )

# Run analysis on sdm data
outcome = "bart"
include("04_analysis.jl")

sdm = (Y = Y,
       richness = richness,
       richness_plot = richness_plot,
       lcbd = lcbd,
       beta_total = beta_total,
       lcbd_plot = lcbdtr_plot,
       rel_plot = rel2d_plot
       )

# Plot distributions
lims_richness = extrema(mapreduce(collect, vcat, [raw.richness, sdm.richness]))
lims_lcbd = extrema(mapreduce(collect, vcat, [raw.lcbd, sdm.lcbd]))

plot(raw.richness_plot, sdm.richness_plot,
     raw.lcbd_plot, sdm.lcbd_plot,
     clim = [lims_richness lims_richness lims_lcbd lims_lcbd],
     title = ["a" "b" "c" "d"],
     titleloc = :left,
     size = (850, 600), 
     dpi = 200)

# Save result
if (@isdefined save_figures) && save_figures == true
    mspath = abspath("..", "ms_betadiversity_hotspots")
    savefig(joinpath(mspath, "figures", "combined-maps.png"))
end

## Correlation, GLM, and friends
# Prepare data
results = DataFrame([raw.richness, raw.lcbd, sdm.richness, sdm.lcbd])
rename!(results, 
        :x1 => :richness_raw, :x2 => :lcbd_raw, 
        :x3 => :richness_sdm, :x4 => :lcbd_sdm)
allowmissing!(results, Not([:latitude, :longitude]))
for col in eachcol(results)
    replace!(col, nothing => missing)
end
dropmissing!(results)
results = convert.(Float64, results)

# Check correlation
cor(results.richness_raw, results.richness_sdm)
cor(results.lcbd_raw, results.lcbd_sdm)

# Check distribution
histogram(results.richness_raw, title = "richness_raw", legend = :none)
histogram(results.richness_sdm, title = "richness_sdm", legend = :none)
histogram(results.lcbd_raw, title = "lcbd_raw", legend = :none)
histogram(results.lcbd_sdm, title = "lcbd_sdm", legend = :none)

mean(results.richness_raw)
std(results.richness_raw)
variation(results.richness_raw)
var(results.richness_raw)/mean(results.richness_raw)

mean(results.richness_sdm)
std(results.richness_sdm)
variation(results.richness_sdm)

## Test linear regression in Julia
using GLM
lm_richness = lm(@formula(richness_sdm ~ richness_raw), results)
lm_lcbd = lm(@formula(lcbd_sdm ~ lcbd_raw), results)

# Test utility functions
coef(lm_richness)
deviance(lm_richness)
dof_residual(lm_richness)
r2(lm_richness)
stderror(lm_richness)
vcov(lm_richness)
predict(lm_richness)

# Test GLMs
glm_richness = glm(@formula(richness_sdm ~ richness_raw), results, Poisson())
negb_richness = negbin(@formula(richness_sdm ~ richness_raw), results, LogLink())
glm_lcbd = glm(@formula(richness_sdm ~ richness_raw), results, Gamma())

## Test regression in R
@rput results

R"""
# Regression
lm_richness <- lm(richness_sdm ~ richness_raw, data = results)
lm_lcbd <- lm(lcbd_sdm ~ lcbd_raw, data = results)

summary(lm_richness)
summary(lm_lcbd)

# Check for assumptions
plotlm <- function(model) {
    opar <- par(mfrow=c(2,2))
    plot(model)
    par(opar)    
}
plotlm(lm_richness) # not met
plotlm(lm_lcbd) # not met

# Plot relation
plot(results$richness_raw, results$richness_sdm)
abline(lm_richness, col = "red")

plot(results$lcbd_raw, results$lcbd_sdm)
abline(lm_lcbd, col = "red")

# Check distribution
hist(results$richness_raw)
hist(results$richness_sdm)
hist(results$lcbd_raw)
hist(results$lcbd_sdm)

# Check transformation
# loglm_richness <- lm(log10(richness_sdm) ~ log10(richness_raw), data = results)
# summary(loglm_richness)
# plotlm(loglm_richness)

# loglm_lcbd <- lm(log10(lcbd_sdm) ~ log10(lcbd_raw), data = results)
# summary(loglm_lcbd)
# plotlm(loglm_lcbd)

# Richness GLMs
# Poisson
glm_richness <- glm(richness_sdm ~ richness_raw, data = results, family = poisson)
summary(glm_richness)
# Overdispersion (Null deviance ~3-4x degrees of freedom)
# Quasi-Poisson
glm_richness <- glm(richness_sdm ~ richness_raw, data = results, family = quasipoisson)
summary(glm_richness)
# Negative binomial
library(MASS)
glm_nb_richness <- glm.nb(richness_sdm ~ richness_raw, data = results)
summary(glm_nb_richness)

# LCBD GLM
glm_lcbd <- glm(lcbd_sdm ~ lcbd_raw, data = results, family = Gamma)
summary(glm_lcbd)

# Get residuals
richness_res <- residuals(glm_nb_richness)
lcbd_res <- residuals(glm_lcbd)

"""

@rget richness_res lcbd_res


## Residual visualization
# Arrange data
results.richness_res = richness_res
results.lcbd_res = lcbd_res
richres_layer = SimpleSDMResponse(results, :richness_res, similar(raw.richness), 
                                  latitude = :latitude, longitude = :longitude)
lcbdres_layer = SimpleSDMResponse(results, :lcbd_res, similar(raw.lcbd), 
                                  latitude = :latitude, longitude = :longitude)

# Plot residuals
plotSDM2(richres_layer, c = :PuOr, dpi = 200)
plotSDM2(lcbdres_layer, c = :PuOr, dpi = 200)
# Check distribution
histogram(richres_layer)
histogram(lcbdres_layer)
histogram2d(richres_layer, lcbdres_layer)

# Julia residuals
histogram(residuals(lm_richness))
histogram(residuals(lm_lcbd))

results.richness_res_jl = residuals(lm_richness)
results.lcbd_res_jl = residuals(lm_lcbd)

richresjl_layer = SimpleSDMResponse(results, :richness_res_jl, similar(raw.richness), 
                                  latitude = :latitude, longitude = :longitude)
lcbdresjl_layer = SimpleSDMResponse(results, :lcbd_res_jl, similar(raw.lcbd), 
                                  latitude = :latitude, longitude = :longitude)

plotSDM2(richresjl_layer, c = :PuOr, dpi = 200)
plotSDM2(lcbdresjl_layer, c = :PuOr, dpi = 200)

