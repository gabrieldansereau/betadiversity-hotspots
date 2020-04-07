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

function random_idx_subset(idx, ratio)
	Random.seed!(42)
	shuffle!(idx)
	ntraining = convert(Int64, round(ratio*size(idx, 1); digits=0))
	trn_idx = idx[1:ntraining]
	vld_idx = idx[ntraining+1:end]
	return (trn = trn_idx, vld = vld_idx)
end

function separate_trn_vld(sp)
	idx0 = findall(iszero, sp)
	idx1 = findall(isone, sp)

	idx0_trn, idx0_vld = random_idx_subset(idx0, 0.7)
	idx1_trn, idx1_vld = random_idx_subset(idx1, 0.7)

	idx_trn = vcat(idx0_trn, idx1_trn)
	idx_vld = vcat(idx0_vld, idx1_vld)
	return (trn = idx_trn, vld = idx_vld)
end

function get_labels(labels)
	trn_idx, vld_idx = separate_trn_vld(labels)
	trn_labels, vld_labels = labels[trn_idx], labels[vld_idx]
	return trn_labels, vld_labels
end

trn_labels, vld_labels = [get_labels(sp) for sp in eachcol(spe)]
