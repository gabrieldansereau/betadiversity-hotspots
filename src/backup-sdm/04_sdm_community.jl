import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load predictions for all species
@load "data/jld2/sdm-predictions-landcover.jld2" predictions

## Create custom functions
# Function for Pielou's evenness index
function pielou(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ sum(A)
    return abs(sum(p.*log.(p))/length(p))
end
function pielou_allspecies(a::Vector{T}) where {T <: Number}
    A = filter(!isnan, a)
    length(A) == 0 && return NaN
    sum(A) == zero(T) && return NaN
    p = A ./ length(A)
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
output2 = zeros(Float64, size(predictions[1]))
# Loop for each pixel/grid element
@time for i in 1:size(output, 1), j in 1:size(output, 2)
    # Group predictions for all species in pixel [i,j]
    x = map(x -> x.grid[i,j], predictions)
    # Calculate Pielou's evenness index for pixel [i,j]
    output[i,j] = pielou(x)
    output2[i,j] = pielou_allspecies(x)
end
# Create SimpleSDMLayer with diversity/evenness scores
diversity = SimpleSDMResponse(output, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)
diversity2 = SimpleSDMResponse(output2, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot result
diversity_plot = plotSDM(diversity, c=:BuPu)
title!(diversity_plot, "SDM species diversity (Pielou's evenness index - Site richness)")
diversity_plot2 = plotSDM(diversity2, c=:BuPu)
title!(diversity_plot2, "SDM species diversity (Pielou's evenness index - Total richness)")

## Save result
#=
savefig(diversity_plot, "fig/sdm/04_sdm_diversity-pielou.pdf")
savefig(diversity_plot2, "fig/sdm/04_sdm_diversity-pielou2.pdf")
=#
