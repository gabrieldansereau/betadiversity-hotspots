library(party)

spe <-  read.csv("data/proc/distributions_spe_qc.csv", header=TRUE, sep="\t")
spa <-  read.csv("data/proc/distributions_spa_qc.csv", header=TRUE, sep="\t")
env <-  read.csv("data/proc/distributions_env_qc.csv", header=TRUE, sep="\t")

head(spe)

set.seed(42)
vars <- cbind(env, spa)
rf <- cforest(as.factor(spe$sp20) ~ .,
              data = vars,
              control = cforest_unbiased(mtry = 2, ntree = 100))

set.seed(42)
(vars_imp <- varimp(rf, conditional = T))
(vars_imp <- varimp(rf))

set.seed(42)
sp1_pred <- predict(rf, OOB=TRUE)

table(spe$sp1, sp1_pred)
