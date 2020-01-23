import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

## Load presence-absence data for all species
@load "data/jld2/raw-pres-abs.jld2" pres_abs

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
    x = map(x -> x.grid[i, j], pres_abs)
    # Calculate Pielou's evenness index for pixel [i,j]
    output[i,j] = pielou(x)
    output2[i,j] = pielou_allspecies(x)
end
# Create SimpleSDMLayer with diversity/evenness scores
diversity = SimpleSDMResponse(output, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)
diversity2 = SimpleSDMResponse(output2, pres_abs[1].left, pres_abs[1].right, pres_abs[1].bottom, pres_abs[1].top)

## Plot result
diversity_plot = plotSDM(diversity, c=:BuPu)
title!(diversity_plot, "Species diversity (Pielou's evenness index - Site richness)")
diversity_plot2 = plotSDM(diversity2, c=:BuPu)
title!(diversity_plot2, "Species diversity (Pielou's evenness index - Total richness)")

## Save result
#=
savefig(diversity_plot, "fig/raw/04_raw_diversity-pielou.pdf")
savefig(diversity_plot2, "fig/raw/04_raw_diversity-pielou2.pdf")
=#
