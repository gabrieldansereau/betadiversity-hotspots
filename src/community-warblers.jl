using Distributed
using JLD2
@time include("required.jl")

## Get & prepare occurence data
@time begin
    # Load data from CSV files
    df = CSV.read("../data/warblers_can.csv", header=true, delim="\t")
    # Prepare data (select columns, arrange values)
    df = prepare_csvdata(df)
    # Separate species
    taxa_occ = [df[df.species .== u,:] for u in unique(df.species)]
end

## Load predictions for all species
@load "../data/predictions-can.jld2" predictions

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
evenness = SDMLayer(output, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot result
sdm_plot = plotSDM(evenness, type="sdm")

## Save result
savefig(sdm_plot, "fig/warblers/warblers-can.pdf")
