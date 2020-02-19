import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load data
spe = CSV.read("data/proc/distributions_spe.csv", header=true, delim="\t")
spa = CSV.read("data/proc/distributions_spa.csv", header=true, delim="\t")
env = CSV.read("data/proc/distributions_env.csv", header=true, delim="\t")

## Perform RandomForests
using RCall
@rput spe spa env
begin
    R"""
    library(party)
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
    """
end
