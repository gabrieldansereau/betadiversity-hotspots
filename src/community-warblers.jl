using Distributed
using JLD2
@time include("required.jl")

## Load predictions for all species
@load "../data/predictions-ebd.jld2" predictions

## Create custom functions
# Function for Pielou's evenness index
function pielou(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ sum(A)
    return abs(sum(p.*log.(p))/length(p))
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
output = zeros(Float64, size(predictions[1]))
# Loop for each pixel/grid element
@time for i in 1:size(output, 1), j in 1:size(output, 2)
    # Group predictions for all species in pixel [i,j]
    x = getindex.(predictions, i, j)
    # Calculate Shannon diversity index for pixel [i,j]
    output[i,j] = shannon(x)
end
# Create SDMLayer with diversity/evenness scores
diversity = SDMLayer(output, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot result
diversity_plot = plotSDM(diversity, type="sdm")
title!(diversity_plot, "Species diversity (Shannon diversity index)")

## Save result
# savefig(diversity_plot, "fig/sdm-diversity.pdf")
