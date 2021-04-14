source(file.path("src", "required.R"))

results <- read_tsv(here("data", "proc", "comparison-results.csv"))
results

lm_richness <- lm(richness_sdm ~ richness_raw, data = results)
lm_lcbd <- lm(lcbd_sdm ~ lcbd_raw, data = results)

summary(lm_richness)
summary(lm_lcbd)

tidy(lm_richness)
glance(lm_richness)
augment(lm_richness)

tidy(lm_lcbd)
glance(lm_lcbd)
augment(lm_lcbd)

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

## Richness GLMs

# Poisson
glm_richness <- glm(richness_sdm ~ richness_raw, data = results, family = poisson)
summary(glm_richness)
tidy(glm_richness)
glance(glm_richness)
# Overdispersion (Null deviance ~3-4x degrees of freedom)

# Quasi-Poisson
glm_qp_richness <- glm(richness_sdm ~ richness_raw, data = results, family = quasipoisson)
summary(glm_qp_richness)
tidy(glm_qp_richness)
glance(glm_qp_richness)

# Negative binomial
library(MASS)
glm_nb_richness <- glm.nb(richness_sdm ~ richness_raw, data = results)
summary(glm_nb_richness)
tidy(glm_nb_richness)
glance(glm_nb_richness)

# LCBD GLM
glm_lcbd <- glm(lcbd_sdm ~ lcbd_raw, data = results, family = Gamma)
summary(glm_lcbd)
tidy(glm_lcbd)
glance(glm_lcbd)

# Get residuals
richness_res <- residuals(glm_nb_richness)
lcbd_res <- residuals(glm_lcbd)