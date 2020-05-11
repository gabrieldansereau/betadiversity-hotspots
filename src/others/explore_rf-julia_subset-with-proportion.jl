import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere begin
    include(joinpath("src", "required.jl"))
    using DecisionTree
    include(joinpath("src", "lib", "model-evaluation.jl"))
end

## Load data
spe = CSV.read(joinpath("data", "proc", "distributions_spe.csv"), header=true, delim="\t")
spa = CSV.read(joinpath("data", "proc", "distributions_spa.csv"), header=true, delim="\t")
env = CSV.read(joinpath("data", "proc", "distributions_env.csv"), header=true, delim="\t")
var = hcat(env, spa)

## Create training and validation sets, with partitioning propotional to class frequency
import MLJ.partition
function partition_by_class(values, fraction; shuffle = true, rng = 42)
    # Get indices for each class
    idx0 = findall(iszero, values)
    idx1 = findall(isone, values)

    # Partition each class into trn, vld
    idx0_trn, idx0_vld = partition(idx0, fraction, shuffle = shuffle, rng = rng)
    idx1_trn, idx1_vld = partition(idx1, fraction, shuffle = shuffle, rng = rng)

    # Combine trn, vld partitions
    idx_trn = vcat(idx0_trn, idx1_trn)
    idx_vld = vcat(idx0_vld, idx1_vld)
    return (trn = idx_trn, vld = idx_vld)
end

# Create empty Arrays
trn_idx = Array{Array{Int64,1},1}()
vld_idx = Array{Array{Int64,1},1}()
# Collect training & validation indices
spe = spe[:, Not(60)] # sp60 has no observations...
@time for sp in eachcol(spe)
    part = partition_by_class(sp, 0.7)
    push!(trn_idx, part.trn)
    push!(vld_idx, part.vld)
end

# Get training and validation sets
function partition_labels_features(labels, features, indices)
    labels_set = Array{Array{Int64,1},1}()
    features_set = Array{Array{Float64,2},1}()
    for i in 1:ncol(labels)
        push!(labels_set, Int64.(labels[indices[i], i]))
        push!(features_set, Array(features[indices[i], :]))
    end
    return labels_set, features_set
end
trn_labels, trn_features = partition_labels_features(spe, var, trn_idx);
vld_labels, vld_features = partition_labels_features(spe, var, vld_idx);

## Train Random Forests for all species
# Set parameters
ntrees = 1000
nsubfeatures = -1
nfolds = 3

# Pair corresponding labels and features
trn_pairs = Pair.(trn_labels, trn_features)
vld_pairs = Pair.(vld_labels, vld_features)

# Train models
@time models = [build_forest(lab, feat, nsubfeatures, ntrees) for (lab, feat) in trn_pairs];
# @save joinpath("..", "julia-rf.jld2") models
# @load joinpath("..", "julia-rf.jld2") models

# Apply models
@time predictions = apply_forest.(models, vld_features)
@time predictions_proba = [apply_forest_proba(m, f, [0,1]) for (m, f) in Pair.(models, vld_features)]

## Evaluate model performance
# Get accuracy measures
@time acc_mes = accuracy_measures.(predictions, vld_labels)
# Get AUC measures
@time auc_stats = auc.(models, vld_features, vld_labels)
# Combine evaluation measures
res = DataFrame(spe = string.("sp", eachindex(acc_mes)),
                freq_abs = Int.(map(sum, eachcol(spe))),
                accuracy = map(x -> x.accuracy, acc_mes),
                sensitivity = map(x -> x.sensitivity, acc_mes),
                specificity = map(x -> x.specificity, acc_mes),
                kappa = map(x -> x.kappa, acc_mes),
                auc = map(x -> x.score, auc_stats),
                plot = map(x -> x.plot, auc_stats)
                )

# Explore measures
auc_plot(1 .- res.specificity, res.sensitivity, NaN)
quantile(res.accuracy)
quantile(filter(!isnan, res.sensitivity))
quantile(filter(!isnan, res.specificity))

# Explore measures for models with > 100 obs only
res100 = filter(x -> x.freq_abs > 100, res)
quantile(res100.accuracy)
quantile(filter(!isnan, res100.sensitivity))
quantile(filter(!isnan, res100.specificity))

# Combine all ROC curves in single plot
p = auc_plot([0], [0], NaN, dpi = 200)
for auc in auc_stats
	plot!(p, auc.FP_rates, auc.sensitivities)
end
p
savefig(p, "fig/random-forests/08_rf_auc.png")
