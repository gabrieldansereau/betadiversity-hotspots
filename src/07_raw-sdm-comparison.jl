import Pkg; Pkg.activate(".")
include("required.jl")

## Run analysis on raw data
outcome = "raw"
include("04_analysis.jl")

# Add non-relative LCBD values
lcbdnr = calculate_lcbd(Y, lcbd; relative = false)

# Assemble results
raw = (Y = Y,
       richness = richness,
       richness_plot = richness_plot,
       lcbd = lcbd,
       lcbdnr = lcbdnr,
       beta_total = beta_total,
       lcbd_plot = lcbdtr_plot,
       rel_plot = rel2d_plot
       )

## Run analysis on sdm data
outcome = "bart"
include("04_analysis.jl")

# Add non-relative LCBD values
lcbdnr = calculate_lcbd(Y, lcbd; relative = false)

# Assemble results
sdm = (Y = Y,
       richness = richness,
       richness_plot = richness_plot,
       lcbd = lcbd,
       lcbdnr = lcbdnr,
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
rescalegrad(:PuOr, extrema(richness_diff); rev = true)
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
results = DataFrame([raw.richness, raw.lcbd, sdm.richness, sdm.lcbd, raw.lcbdnr, sdm.lcbdnr])
rename!(results, 
        :x1 => :richness_raw, :x2 => :lcbd_raw, 
        :x3 => :richness_sdm, :x4 => :lcbd_sdm,
        :x5 => :lcbdnr_raw, :x6 => :lcbdnr_sdm)
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

# LM
lm_richness = lm(@formula(richness_sdm ~ richness_raw), results)
lm_lcbd = lm(@formula(lcbd_sdm ~ lcbd_raw), results)

# Test utility functions
coeftable(lm_richness)
coef(lm_richness)
deviance(lm_richness)
dof_residual(lm_richness)
r2(lm_richness)
stderror(lm_richness)
vcov(lm_richness)
predict(lm_richness)

function glance(model::StatsModels.TableRegressionModel{<:LinearModel})
    DataFrame(
        r_squared = r2(model),
        adj_r_squared = adjr2(model),
        # sigma(model),
        sigma = missing,
        # ftest(model),
        statistic = missing,
        # pvalue(model),
        pvalue = missing,
        df = dof(model),
        loglik = loglikelihood(model),
        AIC = aic(model),
        BIC = bic(model),
        deviance = deviance(model),
        df_residual = dof_residual(model),
        nobs = nobs(model),
    )
end
glance(lm_richness)
glance(lm_lcbd)

# Test GLMs
# Poisson
glm_richness = glm(@formula(richness_sdm ~ richness_raw), results, Poisson())
deviance(glm_richness) / dof_residual(glm_richness) # Overdispersion

function glance(model::StatsModels.TableRegressionModel{<:GeneralizedLinearModel})
    glancedf = DataFrame(
        # null_deviance = nulldeviance(model) # not working
        null_deviance = missing,
        df_null = dof_residual(model),
        loglik = loglikelihood(model),
        AIC = aic(model),
        BIC = bic(model),
        deviance = deviance(model),
        df_residual = nobs(model) - dof(model),
        nobs = nobs(model),
    )
end
glance(glm_richness)

# Negative binomial
negb_richness = negbin(@formula(richness_sdm ~ richness_raw), results, LogLink())
glance(negb_richness)

# LCBD Gamma
glm_lcbd = glm(@formula(lcbd_sdm ~ lcbd_raw), results, Gamma())
glance(glm_lcbd)

# residuals(glm_richness)

## Test regression in R
# Export to CSV
CSV.write(joinpath("data", "proc", "comparison-results.csv"), results, delim = "\t")
# results = CSV.read(joinpath("data", "proc", "comparison-results.csv"), DataFrame)

## Residual visualization
# Load residuals
residuals_df = CSV.read(joinpath("data", "proc", "comparison-residuals.csv"), DataFrame)

# Arrange data
richres_layer = SimpleSDMResponse(residuals_df, :richness, similar(raw.richness), 
                                  latitude = :latitude, longitude = :longitude)
richres_qp_layer = SimpleSDMResponse(residuals_df, :richness_qp, similar(raw.richness), 
                                  latitude = :latitude, longitude = :longitude)
richres_nb_layer = SimpleSDMResponse(residuals_df, :richness_nb, similar(raw.richness), 
                                  latitude = :latitude, longitude = :longitude)
lcbdres_layer = SimpleSDMResponse(residuals_df, :lcbd, similar(raw.lcbd), 
                                  latitude = :latitude, longitude = :longitude)
lcbdres_br_layer = SimpleSDMResponse(residuals_df, :lcbd_br, similar(raw.lcbd), 
                                  latitude = :latitude, longitude = :longitude)

# Plot residuals
plotSDM2(richres_layer, c = :PuOr, dpi = 200)
plotSDM2(richres_qp_layer, c = :PuOr, dpi = 200)
plotSDM2(richres_nb_layer, c = :PuOr, dpi = 200)
plotSDM2(lcbdres_layer, c = :PuOr, dpi = 200)
plotSDM2(lcbdres_br_layer, c = :PuOr, dpi = 200)
# Check distribution
histogram(richres_layer)
histogram(richres_qp_layer)
histogram(richres_nb_layer)
histogram(lcbdres_layer)
histogram(lcbdres_br_layer)
histogram2d(richres_layer, lcbdres_layer)
histogram2d(richres_nb_layer, lcbdres_layer)

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

