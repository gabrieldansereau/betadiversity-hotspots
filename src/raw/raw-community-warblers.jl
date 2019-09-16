using Distributed
using JLD2
@time include("required.jl")

## Load presence-absence data for all species
@load "data/jld2/pres-abs-ebd.jld2" pres_abs

## Create custom functions
# Function for Pielou's evenness index
function pielou(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ sum(A)
    return abs(sum(filter(!isnan, p.*log.(p)))/length(p))
end
function pielou_allspecies(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ length(A)
    return abs(sum(filter(!isnan, p.*log.(p)))/length(p))
end

# Function for Shannon's diversity index
function shannon(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ sum(A)
    return abs(sum(p.*log.(p)))
end

## Calculate diversity/evenness scores
# Empty array for diversity scores
output = zeros(Float64, size(pres_abs[1]))
output2 = zeros(Float64, size(pres_abs[1]))
# Loop for each pixel/grid element
@time for i in 1:size(output, 1), j in 1:size(output, 2)
    # Group presence-absence data for all species in pixel [i,j]
    x = getindex.(pres_abs, i, j)
    # Calculate Pielou's evenness index for pixel [i,j]
    output[i,j] = pielou(x)
    output2[i,j] = pielou_allspecies(x)
end
# Create SDMLayer with diversity/evenness scores
diversity = SDMLayer(output, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)
diversity2 = SDMLayer(output2, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot result
diversity_plot = plotSDM(diversity, type="sdm")
title!(diversity_plot, "Species diversity (Pielou's evenness index - Site richness)")
diversity_plot2 = plotSDM(diversity2, type="sdm")
title!(diversity_plot2, "Species diversity (Pielou's evenness index - Total richness)")

## Save result
savefig(diversity_plot, "fig/raw/raw-diversity-pielou.pdf")
savefig(diversity_plot2, "fig/raw/raw-diversity-pielou2.pdf")
