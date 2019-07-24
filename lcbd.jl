using Plots
using GDAL
using Shapefile
using GBIF
using StatsBase
using Statistics

include("lib/SDMLayer.jl")
include("lib/gdal.jl")
include("lib/worldclim.jl")
include("lib/bioclim.jl")
include("lib/shapefiles.jl")

function gbifdata(sp)
    @info sp
    q = Dict{Any,Any}("limit" => 200, "country" => "CA")
    occ = occurrences(sp, q)
    [next!(occ) for i in 1:2]
    qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates])
    return occ
end

high_taxon = taxon("Anseriformes"; strict=false)
high_taxon_latest_200 = occurrences(high_taxon, Dict("limit"=>200, "country"=>"CA"))
canadian_taxa = unique([p.taxon for p in high_taxon_latest_200])

taxa_occ = gbifdata.(canadian_taxa)
lon_range = (-136.0, -58.0)
lat_range = (40.5, 56.0)

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

@time predictions = [species_bclim(w, wc_vars) for w in taxa_occ]

# Get the LCBD
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

sdm_plot

savefig("lcbd-map.png")
