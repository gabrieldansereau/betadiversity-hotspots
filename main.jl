using Plots
using Statistics
using GDAL
using GBIF

# Core type
include("lib/SDMLayer.jl")

# Read geotiff files
include("lib/gdal.jl")

# Format bioclim data
include("lib/worldclim.jl")

# Utility functions for GBIF + example data
function Base.getindex(p::SDMLayer, r::GBIFRecord)
    return p[r.longitude, r.latitude]
end

sp = taxon("Sphyrapicus varius")
occ = occurrences(sp, Dict{Any,Any}("country" => "CA"))
[next!(occ) for i in 1:20]
function is_in_canada(r::GBIFRecord)
    return r.country == "Canada"
end
qualitycontrol!(occ; filters=[have_ok_coordinates, have_both_coordinates, is_in_canada])

# Get the worldclim data by their layer number
lay1 = worldclim(2)

# Brick the layer by species observations
lats = Float64[]
lons = Float64[]
for r in occ
    push!(lats, r.latitude)
    push!(lons, r.longitude)
end
lay1 = lay1[(minimum(lons)-1.0, maximum(lons)+1.0), (minimum(lats)-1.0, maximum(lats)+1.0)]

heatmap(longitudes(lay1), latitudes(lay1), lay1.grid, aspectratio=1.3, c=:viridis)
scatter!(lons, lats, c=:black, msw=0.0, ms=3, lab="")

obs1 = zeros(Float64, count(occ.show))
i = 1
for r in occ
    obs1[i] = lay1[r]
    global i += 1
end
sort!(obs1)

function find_quantile(itr, p; ql=100)
    qrange = range(0.0; stop=1.0, length=ql)
    q = quantile(itr, qrange; sorted=true)
    return qrange[findmin(abs.(q.-p))[2]]
end

# TODO add getindex(::SDMLayer, ::Int64)
lay1quantiles = zeros(Float64, size(lay1.grid))

for i in eachindex(lay1.grid)
    if isnan(lay1.grid[i])
        local_quantile = NaN
    else
        local_quantile = find_quantile(obs1, lay1.grid[i])
        if local_quantile > 0.5
            local_quantile = 1.0-local_quantile
        end
        local_quantile = 2.0 * local_quantile
    end
    lay1quantiles[i] = local_quantile
end

lay1pred = SDMLayer(lay1quantiles, lay1.left, lay1.right, lay1.bottom, lay1.top)

heatmap(longitudes(lay1pred), latitudes(lay1pred), round.(lay1pred.grid; digits=1), aspectratio=1.3, c=:viridis, clim=(0,1))
scatter!(lons, lats, c=:black, msw=0.0, ms=3, lab="")
