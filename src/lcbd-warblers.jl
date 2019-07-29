using Distributed
addprocs(9)

@time @everywhere include("src/required.jl")

## Get data from CSV files
@time @everywhere begin
    df = CSV.read("../data/warblers_qc_2018.csv", header=true, delim="\t")
    df = prepare_csvdata(df)
    taxa_occ = [df[df.species .== u,:] for u in unique(df.species)]

    lon_range = (-136.0, -58.0)
    lat_range = (40.5, 56.0)
end

@time wc_vars = pmap(x -> worldclim(x)[lon_range, lat_range], 1:19);

@time predictions = pmap(x -> species_bclim(x, wc_vars), taxa_occ);

# Get the LCBD
begin
    Y = zeros(Int64, (prod(size(predictions[1])),length(taxa_occ)))
    @progress for gc in eachindex(predictions[1].grid)
        R = map(x -> x.grid[gc], predictions)
        global Y[gc,:] = .!isnan.(R)
    end

    S = (Y .- mean(Y; dims=1)).^2.0
    SStotal = sum(S)
    BDtotal = SStotal / (size(Y,1)-1)
    SSj = sum(S; dims=1)
    SCBDj = SSj ./ SStotal
    SSi = sum(S; dims=2)
    LCBDi = SSi ./ SStotal

    t_lcbd = zeros(Float64, size(predictions[1]))
    LCBDi = LCBDi./maximum(LCBDi)
    for i in eachindex(t_lcbd)
        t_lcbd[i] = Y[i] > 0 ? LCBDi[i] : NaN
    end

    LCBD = SDMLayer(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

    worldmap = clip(worldshape(50), LCBD)

    sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, size=(900,450), frame=:box)
    xaxis!(sdm_plot, (LCBD.left,LCBD.right), "Longitude")
    yaxis!(sdm_plot, (LCBD.bottom,LCBD.top), "Latitude")

    for p in worldmap
        sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
        plot!(sdm_plot, sh, c=:lightgrey, lab="")
    end

    heatmap!(
    sdm_plot,
    longitudes(LCBD), latitudes(LCBD), LCBD.grid,
    aspectratio=92.60/60.75, c=:viridis,
    clim=(0.0, maximum(filter(!isnan, LCBD.grid)))
    )

    for p in worldmap
        xy = map(x -> (x.x, x.y), p.points)
        plot!(sdm_plot, xy, c=:grey, lab="")
    end

    return sdm_plot

end


savefig("lcbd-map-qc2018.pdf")
