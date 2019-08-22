using Distributed
using JLD2

include("required.jl")

## Get & prepare data
@time begin
    # Load data from CSV files
    df = CSV.read("../data/ebd/ebd_warblers_cut.csv", header=true, delim="\t")
    # Prepare data (select columns, arrange values)
    df = prepare_ebd_data(df)
    # Separate species
    warblers_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-145.0, -50.0)
    lat_range = (20.0, 75.0)
end

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x), 1:19);
temp = wc_vars[1]

## Load predictions
@load "../data/predictions-ebd.jld2" predictions
pred = predictions[1]

## expand
p1 = temp
p2 = pred

# Get minimum coordinates
#=
min_lon = min(p1.left, p2.left)
max_lon = max(p1.right, p2.right)
min_lat = min(p1.bottom, p2.bottom)
max_lat = max(p1.top, p2.top)
=#

# Get coordinates
lons_p1 = longitudes(p1)
lats_p1 = latitudes(p1)
lons_p2 = longitudes(p2)
lats_p2 = latitudes(p2)

# Get position of p2 of p1
m_lon = findmin(abs.(minimum(lons_p2) .- lons_p1))[2]
M_lon = findmin(abs.(maximum(lons_p2) .- lons_p1))[2]
m_lat = findmin(abs.(minimum(lats_p2) .- lats_p1))[2]
M_lat = findmin(abs.(maximum(lats_p2) .- lats_p1))[2]
p1[(m_lat:M_lat), (m_lon:M_lon)]
# Check is size matches
size(p2)

# Create new layer
new1 = fill(NaN, size(p1))
new2 = fill(NaN, size(p1))
new1[(m_lat:M_lat), (m_lon:M_lon)] .= p1[(m_lat:M_lat), (m_lon:M_lon)]
new2[(m_lat:M_lat), (m_lon:M_lon)] .= p2.grid
newlayer1 = SDMLayer(new1, p1.left, p1.right, p1.bottom, p1.top)
newlayer2 = SDMLayer(new2, p1.left, p1.right, p1.bottom, p1.top)

# Plot result
plotp1 = plotSDM(p1)
plotp2 = plotSDM(p2)
plotn1 = plotSDM(newlayer1)
plotn2 = plotSDM(newlayer2)

# Combine heatmaps
sdm_plot = plotSDM(newlayer1)
heatmap!(
    sdm_plot,
    longitudes(newlayer2), latitudes(newlayer2), # layer range
    newlayer2.grid, # evenness values
    aspectratio=92.60/60.75, # aspect ratio
    c=:viridis, # ~color palette
    clim=(0.0, maximum(filter(!isnan, newlayer2.grid))) # colorbar limits
)
plot!(
    sdm_plot,
    xlims=(p2.left, p2.right),
    ylims=(p2.bottom, p2.top),
    aspectratio=92.60/60.75
)
# Compare to p2
plotp2
