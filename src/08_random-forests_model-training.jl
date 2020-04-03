import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load data
spe = CSV.read("data/proc/distributions_spe_qc.csv", header=true, delim="\t")
spa = CSV.read("data/proc/distributions_spa_qc.csv", header=true, delim="\t")
env = CSV.read("data/proc/distributions_env_qc.csv", header=true, delim="\t")
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

predictions = apply_forest(model, vld_features)
predictions_proba = apply_forest_proba(model, vld_features, [0, 1])

function pred_accuracy(preds, obsv)
    sum(preds .== obsv)/size(preds, 1)
end
pred_accuracy(predictions, vld_labels[:, nsp])

countmap(obsv)
tmp = confusion_matrix(vld_labels[:, nsp], predictions)

function accuracy_measures(obsv, preds)
    cm = confusion_matrix(obsv, preds)
    n = sum(cm.matrix)
    TN, FN, FP, TP = cm.matrix
    sensitivity = TP/(TP + FN)
    FN_rate = 1 - sensitivity
    specificity = TN/(TN + FP)
    FP_rate = 1 - specificity
    percent_correct_classification = (TP + TN)/n
    positive_predictive_power = TP/(TP + FP)
    odds_ratio = (TP * TN)/(FP * FN)
    kappa = ( (TP+TN) - (((TP+FN)*(TP+FP) + (FP+TN)*(FN+TN))/n)) / (n-(((TP+FN)*(TP+FP) + (FP+TN)*(FN+TN))/n))
    TSS = 1 - (sensitivity + specificity)
    acc_mes = (cm = cm.matrix,
               TN = TN,
               FN = FN,
               FP = FP,
               TP = TP,
               sensitivity = sensitivity,
               FN_rate = FN_rate,
               specificity = specificity,
               FP_rate = FP_rate,
               percent_correct_classification = percent_correct_classification,
               positive_predictive_power = positive_predictive_power,
               odds_ratio = odds_ratio,
               kappa = kappa,
               TSS = TSS
               )
    return acc_mes
end
@time tmp = accuracy_measures(vld_labels[:, nsp], predictions)
tmp.cm
tmp.sensitivity
tmp.specificity
tmp.percent_correct_classification
tmp.kappa
tmp

p = plot([0,1], [0,1], xlim=(0,1), ylim=(0,1), aspect_ratio=1)
p = scatter!(p, [tmp.FP_rate], [tmp.sensitivity])

function auc_plot(model, vld_features, vld_labels)
    presence_proba = apply_forest_proba(model, vld_features, [0,1])[:,2]
    thresholds = collect(0.0:0.01:1.0)
    thrsh_preds = [replace(prob -> prob .>= thrsh ? 1 : 0, presence_proba) for thrsh in thresholds]
    acc_mes = [accuracy_measures(vld_labels, pred) for pred in thrsh_preds]
    FP_rates = map(x -> x.FP_rate, acc_mes)
    sensitivities = map(x -> x.sensitivity, acc_mes)
    p = plot([0,1], [0,1],
             xlim=(0,1), ylim=(0,1),
             title = "Receiver operating curve",
             xlabel = "False-positive rate (1 - Specificity)",
             ylabel = "True-positive rate (Sensitivity)" ,
             legend = :none, aspect_ratio=1)
    p = plot!(p, FP_rates, sensitivities)
    return p
end
auc_plot(model, vld_features, vld_labels[:, nsp])
