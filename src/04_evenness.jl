import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
# outcome = "raw" # desired outcome (required)
# outcome = "sdm" # desired outcome (required)
# save_figures = true # should figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
  @warn "'outcome' not defined, must be either 'raw' or 'sdm'"
elseif (outcome != "raw" && outcome != "sdm")
  @warn "'outcome' invalid, must be either 'raw' or 'sdm'"
else
  @info "'outcome' currently set to '$(outcome)'"
end

## Load distributions for all species
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions spenames speindex

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
output = zeros(Float64, size(distributions[1]))
output2 = zeros(Float64, size(distributions[1]))
# Loop for each pixel/grid element
@time for i in 1:size(output, 1), j in 1:size(output, 2)
    # Group distributions for all species in pixel [i,j]
    x = map(x -> x.grid[i,j], distributions)
    # Calculate Pielou's evenness index for pixel [i,j]
    output[i,j] = pielou(x)
    output2[i,j] = pielou_allspecies(x)
end
# Create SimpleSDMLayer with diversity/evenness scores
diversity = SimpleSDMResponse(output, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)
diversity2 = SimpleSDMResponse(output2, distributions[1].left, distributions[1].right, distributions[1].bottom, distributions[1].top)

## Plot result
diversity_plot = plotSDM(diversity, c=:BuPu,
                         title = "Community evenness ($(outcome) distributions)",
                         colorbar_title = "Pielou's evenness index (site richness)",
                         dpi = 300)
diversity_plot2 = plotSDM(diversity2, c=:BuPu,
                          title = "Community evenness ($(outcome) distributions)",
                          colorbar_title = "Pielou's evenness index (total richness)",
                          dpi = 300)

## Save result
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) evenness)"
    savefig(diversity_plot, joinpath("fig", outcome, "04_$(outcome)_diversity-pielou.png"))
    savefig(diversity_plot2, joinpath("fig", outcome, "04_$(outcome)_diversity-pielou2.png"))
else
    @info "Figures not saved ($(outcome) evenness)"
end
