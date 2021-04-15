source(file.path("src", "required.R"))

results <- read_tsv(here("data", "proc", "comparison-results.csv"))
results

lm_richness <- lm(richness_sdm ~ richness_raw, data = results)
lm_lcbd <- lm(lcbd_sdm ~ lcbd_raw, data = results)
lm_lcbdnr <- lm(lcbdnr_sdm ~ lcbdnr_raw, data = results)

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
plotfit <- function(model) {
    model %>% 
        augment(type.predict = "response") %>% 
        rename(.y = 1, .x = 2) %>% 
        ggplot(aes(x = .x, y = .y)) +
            geom_point() +
            # geom_abline(
            #     intercept = coef(lm_richness)[1],
            #     slope = coef(lm_richness)[2]
            # )
            geom_line(
                aes(y = .fitted),
                color = "red"
            ) +
            labs(y = names(model$model)[1], x = names(model$model)[2])
}
plotfit(lm_richness)
plotfit(lm_lcbd)
plotfit(lm_lcbdnr)

# Check distribution
hist(results$richness_raw)
hist(results$richness_sdm)
hist(results$lcbd_raw)
hist(results$lcbd_sdm)

## Richness GLMs

# Poisson
glm_richness <- glm(richness_sdm ~ richness_raw, data = results, family = poisson)
summary(glm_richness)
tidy(glm_richness)
glance(glm_richness)
glance(glm_richness) %>% 
    transmute(dispersion = deviance/df.residual)
# Overdispersion (residual deviance ~3-4x residual degrees of freedom)
# Too much dispersion for Poisson
var(results$richness_raw)/mean(results$richness_raw)
var(results$richness_sdm)/mean(results$richness_sdm)
# Variance way bigger than mean, too big for Poisson distribution

plotfit(glm_richness)

# install.packages("ggiraph")
# install.packages("ggiraphExtra")
install.packages("ggeffects")
# install.packages("ggfortify")

library(ggiraph)
library(ggiraphExtra)
ggPredict(glm_richness, se = TRUE)

library(ggeffects)
plot(ggpredict(glm_richness), add.data = TRUE)

library(ggfortify)
autoplot(glm_richness)

results %>% 
    ggplot(aes(longitude, latitude, colour = richness_raw)) +
    scale_color_viridis() +
    geom_point(size = 1)

results %>% 
    ggplot(aes(longitude, latitude, colour = residuals(glm_richness))) +
    scale_color_gradient2() +
    geom_point(size = 1)

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

cv(results$lcbd_raw)
cv(results$lcbd_sdm)
# Gamma supposedly useful when coefficient of variation > 50%,
# which is not the case now
# Maybe Gamma is not right...

AIC(glm_lcbd, lm_lcbd)
# Gamma GLM has lower AIC

plotfit(glm_lcbd)

# Beta regression attempt
library(betareg)
beta_lcbd <- betareg(lcbd_sdm ~ lcbd_raw, data = results, link = "logit")
beta_lcbdnr <- betareg(lcbdnr_sdm ~ lcbdnr_raw, data = results)
summary(beta_lcbd)
tidy(beta_lcbd)
glance(beta_lcbd)
plotlm(beta_lcbd)
plotfit(beta_lcbd)q

lm_uniqueness <- lm(lcbd_raw ~ richness_raw, data = results)
summary(lm_uniqueness)
plotlm(lm_uniqueness)
AIC(lm_uniqueness)

results <- results %>% 
    mutate(
        richness_raw_std = richness_raw/max(richness_raw),
        richness_sdm_std = richness_sdm/max(richness_sdm)
    )

beta_uniqueness <- betareg(lcbd_sdm ~ richness_sdm, data = results)
summary(beta_uniqueness)
plotlm(beta_uniqueness)
AIC(beta_uniqueness)

beta_uniqueness %>% 
    augment(type.predict = "response") %>% 
    rename(.y = 1, .x = 2) %>% 
    ggplot(aes(x = .x, y = .y)) +
        geom_point() +
        geom_line(
            aes(y = .fitted),
            color = "red"
        )

# Get residuals
residuals_df <- tibble(
    longitude = results$longitude,
    latitude = results$latitude,
    richness_res = residuals(glm_richness),
    richness_qp_res = residuals(glm_qp_richness),
    richness_nb_res = residuals(glm_nb_richness),
    lcbd_res = residuals(glm_lcbd),
)
residuals_df

# Export
write_tsv(residuals_df, here("data", "proc", "comparison-residuals.csv"))

