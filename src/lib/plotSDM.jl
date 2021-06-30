## Functions to plot SimpleSDMLayer elements more easily

# Plot layer as a heatmap with worldmap background
"""
    plot_layer(layer::SimpleSDMLayer; shape=true, kw...)

Plots the layer over a grey background based on a WorldClim layer. Note that the
warning about multiple series sharing a color bar is normal. If `shape=true`, it
will also draw a shape countour around the land based on the shapefile return by
`worldmap` (this is slower). Normal plot options can be passed afterwards.
"""
function plot_layer(layer::SimpleSDMLayer; shape=true, kw...)
    ## Arguments
    # layer: SimpleSDMLayer to plot
    # scatter: add observations as points in scatter plot, requires to define occ
    # occ: observations to represent if scatter=true
    # kw: optional plotting arguments

    # Create background layer
    coords = boundingbox(layer)
    baselayer = similar(SimpleSDMPredictor(WorldClim, BioClim, 1; coords...))

    # Create empty plot
    sdm_plot = heatmap(baselayer; c=:lightgrey, lab="", msw=0.0, ms=0.0, frame=:box)

    # Add SDM output as heatmap
    heatmap!(
        sdm_plot,
        layer;
        legend=true,
        clim=extrema(layer), # colorbar limits
        kw..., # additional keyword arguments
    )

    # Add shape countour around land
    if shape
        # Load & clip worldmap background to SimpleSDMLayer (from shp in /assets folder)
        worldmap = clip(worldshape(50), layer)
        # Redraw polygons' outer lines over heatmap values
        for p in worldmap # loop for each polygon
            # Get outer lines coordinates
            xy = map(x -> (x.x, x.y), p.points)
            # Add outer lines to plot
            plot!(sdm_plot, xy; c=:grey, lab="")
        end
    end

    return sdm_plot
end

# Custom plot recipe for layers
@recipe function plot(layer::T) where {T<:SimpleSDMLayer}
    seriestype --> :heatmap
    if get(plotattributes, :seriestype, :heatmap) in [:heatmap, :contour]
        aspect_ratio --> 92.60 / 60.75
        xlims --> (layer.left, layer.right)
        ylims --> (layer.bottom, layer.top)
        xguide --> "Longitude"
        yguide --> "Latitude"
        lg = copy(layer.grid)
        lg[lg .== nothing] .= NaN
        longitudes(layer), latitudes(layer), lg
    elseif get(plotattributes, :seriestype, :histogram) in [:histogram, :density]
        filter(!isnothing, layer.grid)
    end
end
