import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")
@everywhere include("src/lib/model-evaluation.jl")

## Load data
spe = CSV.read("data/proc/distributions_spe.csv", header=true, delim="\t")
spa = CSV.read("data/proc/distributions_spa.csv", header=true, delim="\t")
env = CSV.read("data/proc/distributions_env.csv", header=true, delim="\t")
var = hcat(env, spa)

## Shuffle row indices
idx = collect(1:nrow(var))
Random.seed!(42)
shuffle!(idx)
# Separate training & validation indices
n_training = convert(Int64, round(0.7*size(idx, 1); digits=0))
trn_idx, vld_idx = idx[1:n_training], idx[n_training+1:end]

# Separate training & validation datasets
trn_features, trn_labels = Array(var[trn_idx,:]), Int64.(Array(spe[trn_idx,:]))
vld_features, vld_labels = Array(var[vld_idx,:]), Int64.(Array(spe[vld_idx,:]))

## Train Random Forests model (1 species)
using DecisionTree
nsp = 1
ntrees = 1000
nsubfeatures = -1
nfolds = 3
@time model = build_forest(trn_labels[:, nsp], trn_features, nsubfeatures, ntrees)
@time accuracy = nfoldCV_forest(trn_labels[:, nsp], trn_features, nfolds, nsubfeatures, ntrees)

@time predictions = apply_forest(model, vld_features)
predictions_proba = apply_forest_proba(model, vld_features, [0, 1])

pred_accuracy(predictions, vld_labels[:, nsp])

countmap(vld_labels[:, nsp])
tmp = confusion_matrix(vld_labels[:, nsp], predictions)

@time tmp = accuracy_measures(vld_labels[:, nsp], predictions)
tmp.cm
tmp.sensitivity
tmp.specificity
tmp.accuracy
tmp.kappa
tmp

p = plot([0,1], [0,1], xlim=(0,1), ylim=(0,1), aspect_ratio=1)
p = scatter!(p, [tmp.FP_rate], [tmp.sensitivity])

@time test = auc(model, vld_features, vld_labels[:, nsp])
test.score
test.plot

# cumsum(roc_plot)
# areaplot([0.0, 0.5, 0.1])

## Repeat for all species
@time models = map(x -> build_forest(trn_labels[:,x], trn_features, nsubfeatures, ntrees), 1:size(trn_labels,2));
@time predictions = map(m -> apply_forest(m, vld_features), models)
predictions_mat = hcat(predictions...)

@time acc_mes = map(x -> accuracy_measures(predictions_mat[:,x], vld_labels[:,x]), 1:size(predictions_mat,2));
@time auc_stats = map(x -> auc(models[x], vld_features, vld_labels[:,x]), 1:length(models))
res = DataFrame(spe = string.("sp", eachindex(acc_mes)),
				freq_abs = Int.(map(sum, eachcol(spe))),
				accuracy = map(x -> x.accuracy, acc_mes),
				sensitivity = map(x -> x.sensitivity, acc_mes),
				specificity = map(x -> x.specificity, acc_mes),
				kappa = map(x -> x.kappa, acc_mes),
				auc = map(x -> x.score, auc_stats),
				plot = map(x -> x.plot, auc_stats)
				)

auc_plot(1 .- res.specificity, res.sensitivity, NaN)
quantile(res.accuracy)
quantile(filter(!isnan, res.sensitivity))
quantile(filter(!isnan, res.specificity))
