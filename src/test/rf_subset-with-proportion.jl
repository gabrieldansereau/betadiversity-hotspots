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
spe = spe[:, Not(60)]
@time for sp in eachcol(spe)
	part = partition_by_class(sp, 0.7)
	push!(trn_idx, part.trn)
	push!(vld_idx, part.vld)
end

# Option 3
@time begin
	trn_labels = Array{Array{Int64,1},1}()
	vld_labels = Array{Array{Int64,1},1}()
	trn_features = Array{Array{Float64,2},1}()
	vld_features = Array{Array{Float64,2},1}()
	for i in 1:ncol(spe)
		push!(trn_labels, Int64.(spe[trn_idx[i], i]))
		push!(vld_labels, Int64.(spe[vld_idx[i], i]))
		push!(trn_features, Array(var[trn_idx[i], :]))
		push!(vld_features, Array(var[vld_idx[i], :]))
	end
end;
rm_labs_feats()

# Option 4
function partition_labels_features(labels, features, indices)
	labels_set = Array{Array{Int64,1},1}()
	features_set = Array{Array{Float64,2},1}()
	for i in 1:ncol(labels)
		push!(labels_set, Int64.(labels[indices[i], i]))
		push!(features_set, Array(features[indices[i], :]))
	end
	return labels_set, features_set
end
@time begin
	trn_labels, trn_features = partition_labels_features(spe, var, trn_idx);
	vld_labels, vld_features = partition_labels_features(spe, var, vld_idx);
end;
rm_labs_feats()

function rm_labs_feats()
	global trn_labels = nothing
	global trn_features = nothing
	global vld_labels = nothing
	global vld_features = nothing
end

using DecisionTree
nsp = 1
ntrees = 1000
nsubfeatures = -1
nfolds = 3

trn_pairs = Pair.(trn_labels, trn_features)
vld_pairs = Pair.(vld_labels, vld_features)

@time models = [build_forest(lab, feat, nsubfeatures, ntrees) for (lab, feat) in trn_pairs];
@time predictions = apply_forest.(models, vld_features)
@time predictions_proba = [apply_forest_proba(m, f, [0,1]) for (m, f) in Pair.(models, vld_features)]

@time acc_mes = accuracy_measures.(predictions, vld_labels)
@time auc_stats = auc.(models, vld_features, vld_labels)
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

res100 = filter(x -> x.freq_abs > 100, res)
quantile(res100.accuracy)
quantile(filter(!isnan, res100.sensitivity))
quantile(filter(!isnan, res100.specificity))
