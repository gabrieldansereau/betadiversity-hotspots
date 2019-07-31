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

## Plot results
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

    # Load & clip worldmap background to SDMLayer (from shp in /assets folder)
    worldmap = clip(worldshape(50), evenness)

    # Create empty plot
    sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box)
    # Adjust axes
    xaxis!(sdm_plot, (evenness.left,evenness.right), "Longitude")
    yaxis!(sdm_plot, (evenness.bottom,evenness.top), "Latitude")

    # Add worldmap background
    for p in worldmap # loop for each polygon
        # Construct polygon from points
        sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
        # Add polygon to plot
        plot!(sdm_plot, sh, c=:lightgrey, lab="")
    end

    # Add SDM output as heatmap
    heatmap!(
        sdm_plot,
        longitudes(evenness), latitudes(evenness), # layer range
        evenness.grid, # evenness values
        aspectratio=92.60/60.75, # aspect ratio
        c=:BuPu, # ~color palette
        clim=(0.0, maximum(filter(!isnan, evenness.grid))) # colorbar limits
    )

    # Redraw polygons' outer lines over heatmap values
    for p in worldmap # loop for each polygon
        # Get outer lines coordinates
        xy = map(x -> (x.x, x.y), p.points)
        # Add outer lines to plot
        plot!(sdm_plot, xy, c=:grey, lab="")
    end

    return sdm_plot

end

## Save result
savefig(sdm_plot, "fig/warblers/warblers-qc2018.pdf")
