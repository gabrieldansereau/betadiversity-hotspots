using Plots
using GDAL
using Shapefile
using GBIF
using StatsBase
using Statistics
using EcologicalNetworks
using EcologicalNetworksPlots

include("lib/SDMLayer.jl")
include("lib/gdal.jl")
include("lib/worldclim.jl")
include("lib/bioclim.jl")
include("lib/shapefiles.jl")

function gbifdata(sp)
    @info sp
    q = Dict{Any,Any}("limit" => 200, "country" => "NZ")
    occ = occurrences(sp, q)
    [next!(occ) for i in 1:2]
    qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates])
    return occ
end

N = nz_stream_foodweb()[1]

# Get taxa in GBIF
taxa = Dict{last(eltype(N)),GBIFTaxon}()
@progress for s in species(N)
    @info s
    try
        t = taxon(s; strict=false)
        taxa[s] = t
    catch
        @info "$s has no match"
    end
end

# Get occurrences
occ_data = Dict([taxon => gbifdata(taxon) for taxon in values(taxa)])

filter!(p -> length(p.second) ≥ 10, occ_data)

LAT = vcat(map.(x -> x.latitude, values(occ_data))...)
LON = vcat(map.(x -> x.longitude, values(occ_data))...)

lon_range = (165.0, 180.0)
lat_range = (-30.0, -50.0)

@time wc_vars = [worldclim(i)[lon_range, lat_range] for i in 1:19];

function species_bclim(occ, vars)
    predictions = [bioclim(vars[i], occ) for i in 1:length(vars)];
    prediction = reduce(minimum, predictions);
    for i in eachindex(prediction.grid)
        if prediction.grid[i] == 0.0
            prediction.grid[i] = NaN
        end
    end
    return prediction
end

@time predictions = Dict([occ_d.first => species_bclim(occ_d.second, wc_vars) for occ_d in occ_data])

# Get the LCBD
Y = zeros(Int64, (prod(size(first(predictions).second)), length(predictions)))
@progress for (i, tax) in enumerate(keys(predictions))
    layer = predictions[tax]
    Y[:,i] = vec(.!isnan.(layer.grid))
end

S = (Y .- mean(Y; dims=1)).^2.0
SStotal = sum(S)
BDtotal = SStotal / (size(Y,1)-1)
SSj = sum(S; dims=1)
SCBDj = SSj ./ SStotal
SSi = sum(S; dims=2)
LCBDi = SSi ./ SStotal


# SCBD plotting
scbd_val = Dict{last(eltype(N)),Float64}()
scbd_temp = Dict(zip(keys(predictions), SCBDj))
for s in species(N)
    sp_key = get(taxa, s, NaN)
    if typeof(sp_key) <: GBIFTaxon
        scbd_val[s] = get(scbd_temp, sp_key, 0.0)
    else
        scbd_val[s] = 0.0
    end
end

I0 = initial_random_layout(N)
[graph_layout!(N, I0) for i in 1:2000]
EcologicalNetworksPlots.finish_layout!(I0)

plot(N, I0, nodesize=scbd_val, markercolor=:white, msw=1.0)

savefig("lcbd-network-plot.png")


# LCBD plotting

t_lcbd = zeros(Float64, size(first(predictions).second))
LCBDi = LCBDi./maximum(LCBDi)
for i in eachindex(t_lcbd)
    t_lcbd[i] = Y[i] > 0 ? LCBDi[i] : NaN
end

LCBD = SDMLayer(t_lcbd, first(predictions).second.left, first(predictions).second.right, first(predictions).second.bottom, first(predictions).second.top)

worldmap = clip(worldshape(50), LCBD)

sdm_plot = plot([0.0], lab="", msw=0.0, ms=0.0, frame=:box, dpi=300)
xaxis!(sdm_plot, (LCBD.left,LCBD.right), "Longitude")
yaxis!(sdm_plot, (LCBD.bottom,LCBD.top), "Latitude")

for p in worldmap
    sh = Shape([pp.x for pp in p.points], [pp.y for pp in p.points])
    plot!(sdm_plot, sh, c=:lightgrey, lab="")
end

heatmap!(
    sdm_plot,
    longitudes(LCBD), latitudes(LCBD), LCBD.grid,
    aspectratio=92.60/60.75, c=:magma,
    clim=(0.0, maximum(filter(!isnan, LCBD.grid)))
    )

for p in worldmap
    xy = map(x -> (x.x, x.y), p.points)
    plot!(sdm_plot, xy, c=:grey, lab="")
end

sdm_plot

savefig("lcbd-network-map.png")
