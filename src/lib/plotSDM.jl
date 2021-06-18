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
        clim=(0.0, maximum(filter(!isnothing, layer.grid))); # colorbar limits
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

function plotSDM2(layer::SimpleSDMLayer; kw...)
    ## Arguments
    # layer: SimpleSDMLayer to plot
    # scatter: add observations as points in scatter plot, requires to define occ
    # occ: observations to represent if scatter=true
    # kw: optional plotting arguments

    # Create background layer
    coords = boundingbox(layer)
    baselayer = similar(SimpleSDMPredictor(WorldClim, BioClim, 1; coords...))

    # Create empty plot
    sdm_plot = heatmap(baselayer, c = :lightgrey, lab="", msw=0.0, ms=0.0, frame=:box)

    # Add SDM output as heatmap
    heatmap!(sdm_plot, layer,
             legend = true,
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

function plotlcbd(layer::SimpleSDMLayer, cbtitle::String; kw...)
    p1 = plot(layer, c = :viridis)
    p2 = plot(frame = :none)
    annotate!(p2, 0.5, 0.5, text(cbtitle, 11, :center, 90.0))
    l = @layout [a b{0.01w}]
    plot(p1, p2, layout = l; kw...)
end

@recipe function plot(layer::T) where {T <: SimpleSDMLayer}
    seriestype --> :heatmap
    if get(plotattributes, :seriestype, :heatmap) in [:heatmap, :contour]
        aspect_ratio --> 92.60/60.75
        xlims --> (layer.left, layer.right)
        ylims --> (layer.bottom, layer.top)
        xguide --> "Longitude"
        yguide --> "Latitude"
        lg = copy(layer.grid)
        lg[lg.==nothing] .= NaN
        longitudes(layer), latitudes(layer), lg
    elseif get(plotattributes, :seriestype, :histogram) in [:histogram, :density]
       filter(!isnothing, layer.grid)
    end
 end