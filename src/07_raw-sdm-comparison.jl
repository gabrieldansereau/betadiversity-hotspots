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

## Difference plots
# Difference between values from SDM models & raw observations
richness_diff = sdm.richness - raw.richness
lcbd_diff = sdm.lcbd - raw.lcbd

# Custom fonctions for gradients

rescale(x, m, M) = (x .- minimum(x))./(maximum(x)-minimum(x)).*(M-m).+m

function rescalegrad(grad, lims; kw...) 
    centervalue = abs(lims[1])/(lims[2] - lims[1])
    rescaled_grad = cgrad(grad, centervalue; kw...)
    return rescaled_grad
end

function subsetgrad(grad, lims; kw...)
    lims_range = range(lims..., length = 20) 
    subgrad = cgrad([getindex(cgrad(grad; kw...), lims_range)...])
    return subgrad
end

function recentergrad(grad, lims; kw...)
    lmin = minimum(lims)
    lmax = maximum(lims)
    absmax = max(abs(lmin), lmax)
    if abs(lmin) < lmax
        rlims = rescale([-absmax, lmin, absmax], 0, 1)
        subsetgrad(:PuOr, rlims[2:3]; kw...)
    elseif abs(lmin) > lmax
        rlims = rescale([-absmax, lmax, absmax], 0, 1)
        subsetgrad(:PuOr, rlims[1:2]; kw...)
    end
end

# Test functions
lims = extrema(richness_diff)
rescale([lims[1], 0, lims[2]], 0, 1)
rescalegrad(:PuOr, extrema(layer); rev = true)
subsetgrad(:PuOr, (0.2, 1.0); rev = true)
recentergrad(:PuOr, lims; rev = true)

# Custom function to visualize the difference
function difference_plot(layer::T; title = "") where T <: SimpleSDMLayer
    # Center colorscale at zero instead of midpoint between extremas
    lims = extrema(layer)
    # centervalue = abs(lims[1])/(lims[2] - lims[1])
    # scalevalues = rescale([lims[1], 0.0, lims[2]], 0, 1)
    diff_map = plotSDM2(layer,
                        # c = cgrad(:PuOr, centervalue, rev = true), 
                        # c = rescalegrad(:PuOr, extrema(layer); rev = true),
                        # c = cgrad(:PuOr, scalevalues, rev = true), 
                        c = recentergrad(:PuOr, lims; rev = true),
                        title = "Difference map",
                        colorbar_title = "Difference from observed value")
    diff_hist = histogram([filter(x -> !isnothing(x) && x >= 0, layer.grid), 
                           filter(x -> !isnothing(x) && x < 0, layer.grid)],
                          bins = :rice, c = [:PuOr cgrad(:PuOr; rev = true)], legend = :none,
                          title = "Difference distribution", 
                          xlabel = "Frequency",
                          orientation = :horizontal,
                          )
    diff_title = plot(annotation = (0.5, 0.5, "$(title)"), framestyle = :none)
    l = @layout [t{0.01h}; a{0.6w} b{0.38w}]
    diff_plot = plot(diff_title, diff_map, diff_hist, 
                     size = (850, 400), layout = l,
                     bottommargin = 3.0mm, dpi = 200)
    return diff_plot
end
richness_diffplot = difference_plot(richness_diff; title = "Predicted richness compared to observed richness")
lcbd_diffplot = difference_plot(lcbd_diff; title = "Predicted LCBD compared to observed LCBD")

# Save figures
if (@isdefined save_figures) && save_figures == true
    savefig(richness_diffplot, joinpath("fig", "$(outcome)", "07_$(outcome)_comparison-richness.png"))
    savefig(lcbd_diffplot, joinpath("fig", "$(outcome)", "07_$(outcome)_comparison-lcbd.png"))
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

