using Distributed
addprocs(9)

@time @everywhere include("src/required.jl")

## Get data from CSV files
@time @everywhere begin
    # Load data
    df = CSV.read("../data/warblers_qc.csv", header=true, delim="\t")
    # Prepare data (select columns, arrange values)
    df = prepare_csvdata(df)
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]
    # Select 1 species only
    occ = warblers_occ[1]
end
# with @everywhere: 26.586195 seconds (36.20 M allocations: 1.756 GiB, 5.09% gc time)
# without : 14.836300 seconds (35.36 M allocations: 1.716 GiB, 6.52% gc time)

## Get the worldclim data
@time wc_vars_full = pmap(worldclim, 1:19)

## Create function - SDM & map for 1 species
@everywhere function map_species_distribution(occ; distributed=true, scatter=false)
    # Get the worldclim data by their layer number
    @info "Extract and crop bioclim variables"
    @time wc_vars = pmap(x -> clip(wc_vars_full[x], occ), 1:19, distributed=distributed);
    # Make the prediction for each layer
    @info "Predictions for each layer"
    @time predictions = pmap(x -> bioclim(wc_vars[x], occ), 1:19, distributed=distributed);
    # Make the final prediction by taking the minimum
    @info "Minimum-consensus aggregation"
    @time prediction = reduce(minimum, predictions);
    # Get the threshold for NaN given a percentile
    @info "Threshold estimation"
    @time threshold = first(quantile(prediction[occ], [0.05]))
    @info "5% threshold:\t$(round(threshold; digits=3))"
    # Filter the predictions based on the threshold
    @info "Final prediction filtering"
    @time for i in eachindex(prediction.grid)
        prediction.grid[i] < threshold && (prediction.grid[i] = NaN)
    end

    # Load & clip worldmap background to SDMLayer (from shp in /assets folder)
    # (points with coordinates -> polygons -> background map)
    # is clip() working? --> worldshape(50).shapes is the same, no filter applied...
    worldmap = clip(worldshape(50), prediction)

    # Create empty plot
    sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box,
                    title = first(unique(occ.species)))
    # Adjust axes
    xaxis!(sdm_plot, (prediction.left,prediction.right), "Longitude")
    yaxis!(sdm_plot, (prediction.bottom,prediction.top), "Latitude")

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
        longitudes(prediction), latitudes(prediction), # layer range
        prediction.grid, # SDM values
        aspectratio=92.60/60.75, # aspect ratio
        c=:BuPu, # ~color palette
        clim=(0.0, maximum(filter(!isnan, prediction.grid))) # colorbar limits
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
            c=:black, msw=0.0, ms=0.1, lab="",
            alpha=0.5
            )
    end

    return sdm_plot

end

##  Map 1 species
@time map_species_distribution(warblers_occ[1])
@time [map_species_distribution(occ) for occ in warblers_occ[1:1]]
savefig("fig/warblers/sdm-warbler-qc-$(first(unique(warblers_occ[1].species))).pdf")

## Map all species
@time maps = pmap(x -> map_species_distribution(x, distributed=false), warblers_occ)



######

## Benchmarks
# Option 1: Loop parallelized function
@time maps = [map_species_distribution(occ) for occ in warblers_occ[1:10]]
# [1:2] 1st call: 36.082691 seconds (65.09 M allocations: 4.166 GiB, 6.21% gc time)
# [1:2] 2nd call: 11.754542 seconds (10.38 M allocations: 1.486 GiB, 3.46% gc time)
# [1:10] 1st call: 50.553191 seconds (47.93 M allocations: 5.620 GiB, 4.84% gc time)
# [1:10] 2nd call: 50.521261 seconds (47.71 M allocations: 5.608 GiB, 5.42% gc time)
# 75.086243 seconds (47.89 M allocations: 5.618 GiB, 3.80% gc time)

# Option 2: Parallelize regular function
@time maps = pmap(x -> map_species_distribution(x, distributed=false), warblers_occ[1:10])
# [1:2] 1st call: 46.308213 seconds (15.96 M allocations: 807.545 MiB, 1.24% gc time)
# [1:2] 2nd call: 43.030357 seconds (1.04 M allocations: 69.684 MiB, 0.10% gc time)
# [1:10] 1st call: 48.106559 seconds (4.44 M allocations: 310.167 MiB, 1.65% gc time)
# [1:10] 2nd call: 20.699065 seconds (4.44 M allocations: 310.156 MiB, 3.60% gc time)
# 43.015193 seconds (4.52 M allocations: 314.195 MiB, 2.15% gc time)

# Option 3: Parallelize parallelized function (Overparallelization!)
# NOT WORKING
@time maps = pmap(x -> map_species_distribution(x), warblers_occ[1:2])
# [1:2] 1st call: 59.111134 seconds (1.18 M allocations: 76.521 MiB, 0.11% gc time)
# [1:2] 2nd call: 46.757941 seconds (1.04 M allocations: 69.652 MiB, 0.11% gc time)
# [1:10] NOT WORKING

# Produce maps for all species
@time maps = [map_species_distribution(occ) for occ in warblers_occ]
# 1st call: 199.467498 seconds (239.79 M allocations: 18.215 GiB, 6.15% gc time)
# 2nd call: 174.185456 seconds (178.89 M allocations: 15.235 GiB, 9.83% gc time)
@time maps = pmap(x -> map_species_distribution(x, distributed=false), warblers_occ)
# 1st call: 78.952961 seconds (20.49 M allocations: 1.359 GiB, 3.56% gc time)
# 2nd call: 39.652166 seconds (18.07 M allocations: 1.237 GiB, 6.43% gc time)
