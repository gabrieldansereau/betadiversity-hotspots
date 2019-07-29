using Distributed
addprocs(9)

@time @everywhere include("src/required.jl")

## Get data from CSV files
@time @everywhere begin
    df = CSV.read("../data/warblers_qc.csv", header=true, delim="\t")
    df = prepare_csvdata(df)
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]
    # occ = warblers_occ[1]
end
# with @everywhere: 26.586195 seconds (36.20 M allocations: 1.756 GiB, 5.09% gc time)
# without : 14.836300 seconds (35.36 M allocations: 1.716 GiB, 6.52% gc time)

## Create function
@everywhere function map_species_distribution(occ)
    # Get the worldclim data by their layer number
    @info "Extract and crop bioclim variables"
    @time wc_vars = pmap(x -> clip(worldclim(x), occ), 1:19);
    # Make the prediction for each layer
    @info "Predictions for each layer"
    @time predictions = pmap(x -> bioclim(wc_vars[x], occ), 1:19);
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

    worldmap = clip(worldshape(50), prediction)

    sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box,
                    title = first(unique(occ.species)))
    xaxis!(sdm_plot, (prediction.left,prediction.right), "Longitude")
    yaxis!(sdm_plot, (prediction.bottom,prediction.top), "Latitude")

    for p in worldmap
        sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
        plot!(sdm_plot, sh, c=:lightgrey, lab="")
    end

    heatmap!(
    sdm_plot,
    longitudes(prediction), latitudes(prediction), prediction.grid,
    aspectratio=92.60/60.75, c=:BuPu,
    clim=(0.0, maximum(filter(!isnan, prediction.grid)))
    )

    for p in worldmap
        xy = map(x -> (x.x, x.y), p.points)
        plot!(sdm_plot, xy, c=:grey, lab="")
    end

    sdm_plot

    #=
    scatter!(
    sdm_plot,
    longitudes(occ), latitudes(occ),
    c=:black, msw=0.0, ms=0.1, lab="",
    alpha=0.5
    )
    =#

    return sdm_plot

end

# Test for 1 species
@time maps = pmap(x -> map_species_distribution(x), warblers_occ[1:2])
@time maps = [map_species_distribution(occ) for occ in warblers_occ[1:2]]
maps[1]
maps[2]
# savefig("sdm.png")
savefig("fig/warblers/sdm-warbler-qc-$(first(unique(warblers_occ[1].species))).pdf")
# Produce maps for all species
@time maps = pmap(x -> map_species_distribution(x), warblers_occ)
@time maps = [map_species_distribution(occ) for occ in warblers_occ]
