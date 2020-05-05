function pred_accuracy(preds, obsv)
    sum(preds .== obsv)/size(preds, 1)
end

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

function auc(model::Ensemble{S,T}, vld_features::Array{S}, vld_labels::Array{T,1}; plot = true, kw...) where {S,T}
    presence_proba = apply_forest_proba(model, vld_features, [0,1])[:,2]
    thresholds = collect(0.0:0.01:1.0)
    thrsh_preds = [Int64.(replace(prob -> prob .>= thrsh ? 1 : 0, presence_proba)) for thrsh in thresholds]
    acc_mes = [accuracy_measures(vld_labels, pred) for pred in thrsh_preds]
    FP_rates = map(x -> x.FP_rate, acc_mes)
    sensitivities = map(x -> x.sensitivity, acc_mes)
    score = auc_score(FP_rates, sensitivities)
    if plot
        p = auc_plot(FP_rates, sensitivities, score; kw...)
        return (score = score, FP_rates = FP_rates, sensitivities = sensitivities, plot = p)
    else
        return (score = score, FP_rates = FP_rates, sensitivities = sensitivities)
    end
end
function auc_score(FP_rates, sensitivities)
    FPrev, sensrev = FP_rates[end:-1:1], sensitivities[end:-1:1]
    deltas = FPrev[2:end] .- FPrev[1:end-1]
    areas = deltas .* sensrev[2:end]
    score = sum(areas)
    return score
end
function auc_plot(FP_rates, sensitivities, score; kw...)
    p = plot([0,1], [0,1],
             xlim=(0,1), ylim=(0,1),
             title = "Receiver operating curve",
             xlabel = "False-positive rate (1 - Specificity)",
             ylabel = "True-positive rate (Sensitivity)" ,
             legend = :none, aspect_ratio=1)
plot!(p, FP_rates, sensitivities,
      ann = (0.8, 0.1, "AUC = $(round(score, digits = 3))");
      kw...)
return p
end

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
