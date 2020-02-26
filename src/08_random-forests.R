library(party)

spe <-  read.csv("data/proc/distributions_spe.csv", header=TRUE, sep="\t")
spa <-  read.csv("data/proc/distributions_spa.csv", header=TRUE, sep="\t")
env <-  read.csv("data/proc/distributions_env.csv", header=TRUE, sep="\t")

head(spe)

set.seed(42)
vars <- cbind(env, spa)
rf <- cforest(spe$sp1 ~ .,
              data = vars,
              control = cforest_unbiased(mtry = 2, ntree = 50))

set.seed(42)
vars_imp <- varimp(rf, conditional = T)

set.seed(42)
sp1_pred <- predict(rf, OOB=TRUE)

table(spe$sp1, sp1_pred)
