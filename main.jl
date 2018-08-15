using Plots
using Statistics
using GDAL
using GBIF
using DelimitedFiles

include("lib.jl")

# Extract coordinates
bc_codes = lpad.(1:4, 2, "0")

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

bc_var = bc_vars[2]

quantile_table = zeros(Float64, (tr[2]-bl[2], tr[1]-bl[1]))
@progress for (i, x) in enumerate(bl[2]:1:tr[2]-1)
    for (j, y) in enumerate(bl[1]:1:tr[1]-1)
        val = bc_var[x, y]
        if isnan(val)
            quantile_table[i,j] = NaN
        else
            this_q = find_quantile(predicates[:,1], val)
            this_q = this_q > 0.5 ? 1.0-this_q : this_q
            quantile_table[i,j] = this_q
        end

    end
end

heatmap(2.0 .* quantile_table, c=:viridis)

for (i, ax) in enumerate(x)
    qx = find_quantile(nx, ax)
    qx = qx > 0.5 ? 1-qx : qx
    for (j, ay) in enumerate(y)
        qy = find_quantile(ny, ay)
        qy = qy > 0.5 ? 1-qy : qy
        this_q = minimum([qx, qy])
        S[i,j] = this_q
        S[i,j] = this_q > 0.5 ? 1 - this_q : this_q
    end
end
