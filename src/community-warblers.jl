using Distributed
addprocs(9)

@time @everywhere include("src/required.jl")

## Get data from CSV files
@time @everywhere begin
    df = CSV.read("../data/warblers_qc_2018.csv", header=true, delim="\t")
    df = prepare_csvdata(df)
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    lon_range = (-136.0, -58.0)
    lat_range = (40.5, 56.0)
end

@time wc_vars = pmap(x -> worldclim(x)[lon_range, lat_range], 1:19);

@time predictions = pmap(x -> species_bclim(x, wc_vars), warblers_occ);

begin
    function pielou(a::Vector{T}) where {T <: Number}
        A = filter(!isnan, a)
        length(A) == 0 && return NaN
        sum(A) == zero(T) && return NaN
        p = A ./ sum(A)
        return abs(sum(p.*log.(p))/length(p))
    end

    function shannon(a::Vector{T}) where {T <: Number}
        A = filter(!isnan, a)
        length(A) == 0 && return NaN
        sum(A) == zero(T) && return NaN
        p = A ./ sum(A)
        return abs(sum(p.*log.(p)))
    end

    output = zeros(Float64, size(predictions[1]))
    @time for i in 1:size(output, 1), j in 1:size(output, 2)
        x = getindex.(predictions, i, j)
        output[i,j] = shannon(x)
    end

    evenness = SDMLayer(output, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

    worldmap = clip(worldshape(50), evenness)

    sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box)
    xaxis!(sdm_plot, (evenness.left,evenness.right), "Longitude")
    yaxis!(sdm_plot, (evenness.bottom,evenness.top), "Latitude")

    for p in worldmap
        sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
        plot!(sdm_plot, sh, c=:lightgrey, lab="")
    end

    heatmap!(
    sdm_plot,
    longitudes(evenness), latitudes(evenness), evenness.grid,
    aspectratio=92.60/60.75, c=:BuPu,
    clim=(0.0, maximum(filter(!isnan, evenness.grid)))
    )

    for p in worldmap
        xy = map(x -> (x.x, x.y), p.points)
        plot!(sdm_plot, xy, c=:grey, lab="")
    end

    return sdm_plot
    
end

savefig("warblers-qc2018.pdf")
