import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load data
spe = CSV.read("data/proc/distributions_spe_qc.csv", header=true, delim="\t")
spa = CSV.read("data/proc/distributions_spa_qc.csv", header=true, delim="\t")
env = CSV.read("data/proc/distributions_env_qc.csv", header=true, delim="\t")
var = hcat(env, spa)

function get_features_and_labels(features, labels)
end
features = spe
labels = var

idx = collect(1:nrow(var))
@time begin
	idx0 = map(x -> findall(iszero, x), eachcol(features))
	idx1 = map(x -> findall(isone, x), eachcol(features))

	Random.seed!(42)
	shuffle!.(idx0)
	shuffle!.(idx1)

	n_training0 = [convert(Int64, round(0.7*size(i, 1); digits=0)) for i in idx0]
	n_training1 = [convert(Int64, round(0.7*size(i, 1); digits=0)) for i in idx1]
	trn_idx0 = [i[1:n_training0] for i in idx0]
	trn_idx1 = [i[1:n_training1] for i in idx1]
end

function separate_indices(idx, ratio)
	Random.seed!(42)
	shuffle!(idx)
	ntraining = convert(Int64, round(ratio*size(idx, 1); digits=0))
	trn_idx = idx[1:ntraining]
	vld_idx = idx[ntraining+1:end]
	return (trn = trn_idx, vld = vld_idx)
end

function separate_features_labels_indices(sp)
	idx0 = findall(iszero, sp)
	idx1 = findall(isone, sp)

	idx0_trn, idx0_vld = separate_indices(idx0, 0.7)
	idx1_trn, idx1_vld = separate_indices(idx1, 0.7)

	idx_trn = vcat(idx0_trn, idx1_trn)
	idx_vld = vcat(idx0_vld, idx1_vld)
	return (trn = idx_trn, vld = idx_vld)
end

tmp = map(separate_features_labels_indices, eachcol(spe))

trn_set = [spe[:,i][tmp[i].trn] for i in 1:ncol(spe)]
vld_set = [spe[:,i][tmp[i].vld] for i in 1:ncol(spe)]

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

function pred_accuracy(preds, obsv)
    sum(preds .== obsv)/size(preds, 1)
end
pred_accuracy(predictions, vld_labels[:, nsp])

countmap(vld_labels[:, nsp])
tmp = confusion_matrix(vld_labels[:, nsp], predictions)

function accuracy_measures(obsv::Array{T,1}, preds::Array{T,1}) where {T}
    cm = confusion_matrix(obsv, preds)
    n = sum(cm.matrix)
	if isone(length(cm.classes))
		@warn "Only one class in confusion matrix"
		TN, FN, FP, TP = 0., 0., 0., 1.
	else
    	TN, FN, FP, TP = cm.matrix
	end
    sensitivity = TP/(TP + FN)
    FN_rate = 1 - sensitivity
    specificity = TN/(TN + FP)
    FP_rate = 1 - specificity
    accuracy = (TP + TN)/n
    positive_predictive_power = TP/(TP + FP)
    odds_ratio = (TP * TN)/(FP * FN)
    kappa = ( (TP+TN) - (((TP+FN)*(TP+FP) + (FP+TN)*(FN+TN))/n)) / (n-(((TP+FN)*(TP+FP) + (FP+TN)*(FN+TN))/n))
    TSS = 1 - (sensitivity + specificity)
    acc_mes = (cm = cm.matrix,
               TN = TN, FN = FN, FP = FP, TP = TP,
               sensitivity = sensitivity,
               FN_rate = FN_rate,
               specificity = specificity,
               FP_rate = FP_rate,
               accuracy = accuracy,
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
tmp.accuracy
tmp.kappa
tmp

p = plot([0,1], [0,1], xlim=(0,1), ylim=(0,1), aspect_ratio=1)
p = scatter!(p, [tmp.FP_rate], [tmp.sensitivity])

function auc(model::Ensemble{S,T}, vld_features::Array{S}, vld_labels::Array{T,1}) where {S,T}
	presence_proba = apply_forest_proba(model, vld_features, [0,1])[:,2]
	thresholds = collect(0.0:0.01:1.0)
	thrsh_preds = [Int64.(replace(prob -> prob .>= thrsh ? 1 : 0, presence_proba)) for thrsh in thresholds]
	acc_mes = [accuracy_measures(vld_labels, pred) for pred in thrsh_preds]
	FP_rates = map(x -> x.FP_rate, acc_mes)
	sensitivities = map(x -> x.sensitivity, acc_mes)
	score = auc_score(FP_rates, sensitivities)
 	p = auc_plot(FP_rates, sensitivities, score)
	return (score = score, plot = p)
end
function auc_score(FP_rates, sensitivities)
	FPrev, sensrev = FP_rates[end:-1:1], sensitivities[end:-1:1]
	deltas = FPrev[2:end] .- FPrev[1:end-1]
	areas = deltas .* sensrev[2:end]
	score = sum(areas)
	return score
end
function auc_plot(FP_rates, sensitivities, score)
	p = plot([0,1], [0,1],
			 xlim=(0,1), ylim=(0,1),
			 title = "Receiver operating curve",
			 xlabel = "False-positive rate (1 - Specificity)",
			 ylabel = "True-positive rate (Sensitivity)" ,
			 legend = :none, aspect_ratio=1)
	scatter!(p, FP_rates, sensitivities,
			 ann = (0.8, 0.1, "AUC = $(round(score, digits = 3))"))
	return p
end

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

# From https://github.com/marubontan/MLUtils.jl/blob/master/src/utils.jl
#=
function auc2(yTruth::Array{Int}, yScore::Array{Float64})
	dIndex = findall(1 .== yTruth)
	dnIndex = findall(0 .== yTruth)
	score = 0.0
	for ydn in dnIndex
		for yd in dIndex
			if yScore[yd] > yScore[ydn]
				score += 1
			elseif yScore[ydn] == yScore[yd]
				score += 0.5
			end
		end
	end
	return score / (length(dIndex) * length(dnIndex))
end
yTruth = vld_labels[:, nsp]
yScore = apply_forest_proba(model, vld_features, [0,1])[:,2]
auc2(yTruth, yScore)
=#
