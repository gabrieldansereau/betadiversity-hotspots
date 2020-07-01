## Functions to plot SimpleSDMLayer elements more easily

# Plot layer as a heatmap with worldmap background
function plotSDM(layer::SimpleSDMLayer; scatter::Bool=false, occ=nothing, kw...)
    ## Arguments
    # layer: SimpleSDMLayer to plot
    # scatter: add observations as points in scatter plot, requires to define occ
    # occ: observations to represent if scatter=true
    # kw: optional plotting arguments

    # Load & clip worldmap background to SimpleSDMLayer (from shp in /assets folder)
    worldmap = clip(worldshape(50), layer)

    # Create empty plot
    sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, frame=:box)
    # Adjust axes to layer coordinates
    xaxis!(sdm_plot, (layer.left,layer.right), "Longitude")
    yaxis!(sdm_plot, (layer.bottom,layer.top), "Latitude")

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
        longitudes(layer), latitudes(layer), # layer range
        layer.grid, # layer values
        aspectratio=92.60/60.75, # aspect ratio
        clim=(0.0, maximum(filter(!isnan, layer.grid))); # colorbar limits
        kw... # additional keyword arguments
    )

    # Redraw polygons' outer lines over heatmap values
    for p in worldmap # loop for each polygon
        # Get outer lines coordinates
        xy = map(x -> (x.x, x.y), p.points)
        # Add outer lines to plot
        plot!(sdm_plot, xy, c=:grey, lab="")
    end

    # Add observations as scatter points (if scatter=true in function call)
    if scatter == true
        scatter!(
            sdm_plot,
            longitudes(occ), latitudes(occ),
            c=:black, msw=0.0, ms=2.0,
            ma=0.1, mc=:black,
            lab="",
            alpha=0.5
            )
    end

    return sdm_plot
end
plotSDM(layer::SimpleSDMLayer) = print("LOL XD")

function plotSDM2(layer::SimpleSDMLayer; kw...)
    ## Arguments
    # layer: SimpleSDMLayer to plot
    # scatter: add observations as points in scatter plot, requires to define occ
    # occ: observations to represent if scatter=true
    # kw: optional plotting arguments

    # Create background layer
    baselayer = worldclim(1)[(left = layer.left, right = layer.right,
                              top = layer.top, bottom = layer.bottom)]

    # Create empty plot
    sdm_plot = heatmap(baselayer, c = :lightgrey, lab="", msw=0.0, ms=0.0, frame=:box)

    # Add SDM output as heatmap
    heatmap!(
        sdm_plot,
        longitudes(layer), latitudes(layer), # layer range
        layer.grid, # layer values
        aspectratio = 92.60/60.75, # aspect ratio
        clim = extrema(layer); # colorbar limits
        kw... # additional keyword arguments
    )

    # Load & clip worldmap background to SimpleSDMLayer (from shp in /assets folder)
    worldmap = clip(worldshape(50), layer)
    # Redraw polygons' outer lines over heatmap values
    for p in worldmap # loop for each polygon
        # Get outer lines coordinates
        xy = map(x -> (x.x, x.y), p.points)
        # Add outer lines to plot
        plot!(sdm_plot, xy, c=:grey, lab="")
    end

    return sdm_plot
end

@recipe function plot(layer::T) where {T <: SimpleSDMLayer}
    seriestype --> :heatmap
    @assert eltype(layer) <: Number
    if get(plotattributes, :seriestype, :heatmap) in [:heatmap, :contour]
        aspect_ratio --> 92.60/60.75
        xlims --> (minimum(longitudes(layer)),maximum(longitudes(layer)))
        ylims --> (minimum(latitudes(layer)),maximum(latitudes(layer)))
        xlabel --> "Longitude"
        ylabel --> "Latitude"
        longitudes(layer), latitudes(layer), layer.grid
    elseif get(plotattributes, :seriestype, :histogram) in [:histogram, :density]
        filter(!isnan, layer.grid)
    end
end
