include("required.jl")

## Get data
@load joinpath("data", "jld2", "comparison-results.jld2") raw sdm
results = CSV.read(joinpath("data", "proc", "comparison-results.csv"), DataFrame)

# Check correlation
cor(results.richness_raw, results.richness_sdm)
cor(results.lcbd_raw, results.lcbd_sdm)

# Check distribution
histogram(results.richness_raw; title="richness_raw", legend=:none)
histogram(results.richness_sdm; title="richness_sdm", legend=:none)
histogram(results.lcbd_raw; title="lcbd_raw", legend=:none)
histogram(results.lcbd_sdm; title="lcbd_sdm", legend=:none)

mean(results.richness_raw)
std(results.richness_raw)
variation(results.richness_raw)
var(results.richness_raw) / mean(results.richness_raw)

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
    return DataFrame(;
        r_squared=r2(model),
        adj_r_squared=adjr2(model),
        # sigma(model),
        sigma=missing,
        # ftest(model),
        statistic=missing,
        # pvalue(model),
        pvalue=missing,
        df=dof(model),
        loglik=loglikelihood(model),
        AIC=aic(model),
        BIC=bic(model),
        deviance=deviance(model),
        df_residual=dof_residual(model),
        nobs=nobs(model),
    )
end
glance(lm_richness)
glance(lm_lcbd)

# Test GLMs
# Poisson
glm_richness = glm(@formula(richness_sdm ~ richness_raw), results, Poisson())
deviance(glm_richness) / dof_residual(glm_richness) # Overdispersion

function glance(model::StatsModels.TableRegressionModel{<:GeneralizedLinearModel})
    return glancedf = DataFrame(;
        # null_deviance = nulldeviance(model) # not working
        null_deviance=missing,
        df_null=dof_residual(model),
        loglik=loglikelihood(model),
        AIC=aic(model),
        BIC=bic(model),
        deviance=deviance(model),
        df_residual=nobs(model) - dof(model),
        nobs=nobs(model),
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
