function plotSDM(layer::SDMLayer; type::String="sdm", scatter::Bool=false, occ=nothing)
    ## Arguments
    # layer: SDMLayer to plot
    # type: type of layer to represent, either "sdm" (default) or "lcbd"
    # scatter: add observations as points in scatter plot, requires to define occ
    # occ: observations to represent if scatter=true

    # Load & clip worldmap background to SDMLayer (from shp in /assets folder)
    worldmap = clip(worldshape(50), layer)

    # Create empty plot
    sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box)
    # Adjust axes
    xaxis!(sdm_plot, (layer.left,layer.right), "Longitude")
    yaxis!(sdm_plot, (layer.bottom,layer.top), "Latitude")

    # Add worldmap background
    for p in worldmap # loop for each polygon
        # Construct polygon from points
        sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
        # Add polygon to plot
        plot!(sdm_plot, sh, c=:lightgrey, lab="")
    end

    # Define plot aspect according to type
    if type == "sdm"
        colorpalette = :BuPu
    elseif type == "lcbd"
        colorpalette = :viridis
    end

    # Add SDM output as heatmap
    heatmap!(
        sdm_plot,
        longitudes(layer), latitudes(layer), # layer range
        layer.grid, # evenness values
        aspectratio=92.60/60.75, # aspect ratio
        c=colorpalette, # ~color palette
        clim=(0.0, maximum(filter(!isnan, layer.grid))) # colorbar limits
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