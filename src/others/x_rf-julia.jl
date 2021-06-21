using Distributed
@time begin
    @time include("required.jl")
    using DecisionTree
    include(joinpath("src", "lib", "model-evaluation.jl"))
end

## Load data
spe = CSV.read(
    joinpath("data", "proc", "distributions_spe_full.csv"); header=true, delim="\t"
)
spa = CSV.read(
    joinpath("data", "proc", "distributions_spa_full.csv"); header=true, delim="\t"
)
env = CSV.read(
    joinpath("data", "proc", "distributions_env_full.csv"); header=true, delim="\t"
)

# Select observed sites only
spa = spa[spe.site, :]
env = env[spe.site, :]
var = join(env, spa; on=:site)

# Remove site column
select!(spe, Not(:site))
select!(spa, Not(:site))
select!(env, Not(:site))
select!(var, Not(:site))

## Shuffle row indices
idx = collect(1:nrow(var))
Random.seed!(42)
shuffle!(idx)
# Separate training & validation indices
n_training = convert(Int64, round(0.7 * size(idx, 1); digits=0))
trn_idx, vld_idx = idx[1:n_training], idx[(n_training + 1):end]

# Separate training & validation datasets
trn_features, trn_labels = Array(var[trn_idx, :]), Int64.(Array(spe[trn_idx, :]))
vld_features, vld_labels = Array(var[vld_idx, :]), Int64.(Array(spe[vld_idx, :]))

## Train Random Forests model (1 species)
# Set parameters
nsp = 1
ntrees = 1000
nsubfeatures = -1
nfolds = 3
# Train model
@time model = build_forest(trn_labels[:, nsp], trn_features, nsubfeatures, ntrees)
# Cross validation
@time accuracy = nfoldCV_forest(
    trn_labels[:, nsp], trn_features, nfolds, nsubfeatures, ntrees
)

# Apply model
@time predictions = apply_forest(model, vld_features)
@time predictions_proba = apply_forest_proba(model, vld_features, [0, 1])

# Check prediction accuracy
pred_accuracy(predictions, vld_labels[:, nsp])
# Count number of obs in each class
countmap(vld_labels[:, nsp])
# Get confusion matrix
confusion_matrix(vld_labels[:, nsp], predictions)
# Get accuracy measures
@time acc_mes = accuracy_measures(vld_labels[:, nsp], predictions)
acc_mes.cm
acc_mes.sensitivity
acc_mes.specificity
acc_mes.accuracy
acc_mes.kappa
acc_mes

# Get AUC measures
@time auc_stats = auc(model, vld_features, vld_labels[:, nsp])
auc_stats.score
auc_stats.plot

## Repeat for all species
# Train models
@time models = map(
    x -> build_forest(trn_labels[:, x], trn_features, nsubfeatures, ntrees),
    1:size(trn_labels, 2),
);
# Apply
@time predictions = map(m -> apply_forest(m, vld_features), models)
predictions_mat = hcat(predictions...)

# Get evaluation measures
@time acc_mes = map(
    x -> accuracy_measures(predictions_mat[:, x], vld_labels[:, x]),
    1:size(predictions_mat, 2),
);
@time auc_stats = map(x -> auc(models[x], vld_features, vld_labels[:, x]), 1:length(models))
# Combine measures
res = DataFrame(;
    spe=string.("sp", eachindex(acc_mes)),
    freq_abs=Int.(map(sum, eachcol(spe))),
    accuracy=map(x -> x.accuracy, acc_mes),
    sensitivity=map(x -> x.sensitivity, acc_mes),
    specificity=map(x -> x.specificity, acc_mes),
    kappa=map(x -> x.kappa, acc_mes),
    auc=map(x -> x.score, auc_stats),
    plot=map(x -> x.plot, auc_stats),
)

# Explore measures
auc_plot(1 .- res.specificity, res.sensitivity, NaN)
quantile(res.accuracy)
quantile(filter(!isnan, res.sensitivity))
quantile(filter(!isnan, res.specificity))
