function partition_by_class(sp)
	idx0 = findall(iszero, sp)
	idx1 = findall(isone, sp)

	idx0_trn, idx0_vld = partition(idx0, 0.7, shuffle = true, rng = 42)
	idx1_trn, idx1_vld = partition(idx1, 0.7, shuffle = true, rng = 42)

	idx_trn = vcat(idx0_trn, idx1_trn)
	idx_vld = vcat(idx0_vld, idx1_vld)
	return (trn = idx_trn, vld = idx_vld)
end

function get_labels(labels)
	trn_idx, vld_idx = partition_by_class(labels)
	trn_labels, vld_labels = labels[trn_idx], labels[vld_idx]
	return trn_labels, vld_labels
end

import MLJ.partition
# Create empty Arrays
trn_idx = Array{Array{Int64,1},1}()
vld_idx = Array{Array{Int64,1},1}()
# Collect training & validation indices
@time for sp in eachcol(spe)
	push!(trn_idx, separate_trn_vld(sp).trn)
	push!(vld_idx, separate_trn_vld(sp).vld)
end

@time begin
	# Create labels sets
	trn_labels = [Int64.(spe[trn_idx[i], i]) for i in 1:ncol(spe)];
	vld_labels = [Int64.(spe[vld_idx[i], i]) for i in 1:ncol(spe)];
	# Create features sets
	trn_features = [Array(var[trn_idx[i], :]) for i in 1:ncol(spe)];
	vld_features = [Array(var[trn_idx[i], :]) for i in 1:ncol(spe)];
end;
@time begin
	trn_labels, vld_labels = [map(i -> Int64.(spe[idx[i], i]), 1:ncol(spe)) for idx in (trn_idx, vld_idx)];
	trn_features, vld_features = [map(i -> Array(var[i, :]), idx) for idx in (trn_idx, vld_idx)];
end;
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
end
