using Distributed
addprocs(9)

@time @everywhere include("src/required.jl")

## Get & prepare data
@time @everywhere begin
    # Load data from CSV files
    df = CSV.read("../data/warblers_qc_2018.csv", header=true, delim="\t")
    # Prepare data (select columns, arrange values)
    df = prepare_csvdata(df)
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-136.0, -58.0)
    lat_range = (40.5, 56.0)
end

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x)[lon_range, lat_range], 1:19);

## Make predictions for all species
@time predictions = pmap(x -> species_bclim(x, wc_vars), warblers_occ);

## Calculate diversity/evenness scores
begin
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

    # Calculate diversity index
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
end

## Plot result
sdm_plot = plotSDM(evenness, type="sdm")

## Save result
savefig(sdm_plot, "fig/warblers/warblers-qc2018.pdf")
