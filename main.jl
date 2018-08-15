using Plots
using Statistics
using GDAL
using GBIF
using DelimitedFiles

include("lib.jl")

# Extract coordinates
bc_codes = lpad.(1:19, 2, "0")

bc_vars = [get_tiff_data("assets/wc2.0_bio_10m_$(x).tif") for x in bc_codes]

lat = collect(range(-90.0; stop=90.0, length=size(bc_vars[1], 1)))
lon = collect(range(-180.0; stop=180.0, length=size(bc_vars[1], 2)))
heatmap(lon, lat, bc_vars[4], aspectratio=1, c=:viridis)

# Read occurrences
finch = readdlm("assets/Haemorhous_purpureus.coordinates")
n_obs = size(finch, 1)

predicates = zeros(Float64, (n_obs, length(bc_vars)))
@progress for i in 1:n_obs
    obs_point = (finch[i,1], finch[i,2])
    grid_point = get_closest_grid_point(obs_point, lon, lat)
    for j in 1:length(bc_vars)
        predicates[i,j] = bc_vars[j][reverse(grid_point)...]
    end
end


heatmap(predicates)

heatmap(lon, lat, bc_vars[4])
scatter!(finch[:,1], finch[:,2])
xaxis!((-150,-50))
yaxis!((25,75))

bl = get_closest_grid_point((-130.0,25.0), lon, lat)
tr = get_closest_grid_point((-50.0, 60.0), lon, lat)

function quantile_matrix(bc_var::Matrix{Float64}, obs::Vector{Float64}, lower_left::NTuple{2,Int64}, upper_right::NTuple{2,Int64})
    qtable = zeros(Float64, (upper_right[2]-lower_left[2], upper_right[1]-lower_left[1]))
    @progress for (i, x) in enumerate(lower_left[2]:1:upper_right[2]-1)
        for (j, y) in enumerate(lower_left[1]:1:upper_right[1]-1)
            val = bc_var[x, y]
            if isnan(val)
                qtable[i,j] = NaN
            else
                this_q = find_quantile(obs, val)
                this_q = this_q > 0.5 ? 1.0-this_q : this_q
                qtable[i,j] = this_q
            end

        end
    end
    return 2.0 .* qtable
end

quantiles_by_var = [quantile_matrix(bc_vars[i], predicates[:,i], bl, tr) for i in 1:length(bc_vars)]

consensus_matrix = zeros(Float64, size(quantiles_by_var[1]))
for i in eachindex(consensus_matrix)
    consensus_matrix[i] = minimum(getindex.(quantiles_by_var, i))
end

bbox_lat = lat[bl[2]:tr[2]-1]
bbox_lon = lon[bl[1]:tr[1]-1]
heatmap(bbox_lon, bbox_lat, consensus_matrix, frame=:none, c=:Blues)
scatter!(finch[:,1], finch[:,2], leg=false, msw=0, c=:black, m=:diamond, ms=2)
xaxis!((minimum(bbox_lon), maximum(bbox_lon)))
yaxis!((minimum(bbox_lat), maximum(bbox_lat)))
