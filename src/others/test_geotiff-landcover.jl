import Pkg; Pkg.activate(".")
include("../required.jl")

## Test landcover variables
# Define coordinates range
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)

# Test loading variables
lc_vars = map(x -> landcover(x, resolution = 10.0)[coords], 1:10)

# Temporary fix for landcover layers dimensions
lc_vars = [l[top = coords.top - stride(l, 2)] for l in lc_vars]

# Load worldclim variables to compare
wc_vars = SimpleSDMPredictor(WorldClim, BioClim, 1:19; resolution = 10.0, coords...);

## Load with geotiff instead

# Check which layer is crops
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame) # lc2

# Load crops with geotiff
lc2_geo_full = geotiff(SimpleSDMPredictor, joinpath("assets", "landcover", "lc_crops_10m.tif"))
lc2_geo = geotiff(SimpleSDMPredictor, joinpath("assets", "landcover", "lc_crops_10m.tif"); coords...)
lc2_geo = convert(Float32, lc2_geo)
# Get corresponding layer
lc2_csv = lc_vars[2]

# Tests
isequal(lc2_geo, lc2_csv)
isequal(lc2_geo.grid, lc2_csv.grid)
isapprox(lc2_geo.grid, lc2_csv.grid)
isequal(collect(lc2_geo), collect(lc2_csv))
isequal(length(lc2_geo), length(lc2_csv)) # not same number of elements...

testdf = DataFrame([lc2_geo, lc2_csv]) # not same top coordinates

lc2_geo
lc2_csv
isapprox(stride(lc2_geo, 1), stride(lc2_csv, 1)) # not same stride...

size(lc2_geo_full) # 840 x 2159, but should be 840 x 2160...

## Check details from landcover function

# Arguments
layers = collect(1:10)
resolution = 10.0
path = joinpath("assets", "landcover/")

# Function content
# function landcover(layers::Vector{Int64}; resolution::AbstractFloat=10.0, path::AbstractString=joinpath("assets", "landcover/"))
    ## Get file paths
    # List files in path
    lc_files = readdir("$(path)")
    # Filter for selected resolution
    filter!(x -> occursin.("$(Int(resolution))m.tif", x), lc_files)
    # Create paths for selected layers only
    paths = joinpath.(path, lc_files)[layers]

    ## Load data
    # Read raster dataset
    datasets = [ArchGDAL.read(p) for p in paths]
    datasets[2] # 2160 x 240
    # Read raster values from band 1 (only band in this case)
    values = [ArchGDAL.read(d, 1) for d in datasets]
    values[2] # 2160 x 840
    # HOW???
    # Convert from UInt8 to Float64
    values_int = [convert.(Float64, v) for v in values]

    ## Preliminary manipulations
    # Reverse rows and permute dimensions
    landcover_mat = [permutedims(v[:,end:-1:1]) for v in values_int]
    # Replace 255 (default no data values) by nothing
    landcover_mat = convert.(Array{Union{Float32, Nothing}}, landcover_mat)
    [replace!(l, 255 => nothing) for l in landcover_mat]

    # Fill missing latitudes with nothing (latitude extent is only (-60,80) for landcover data)
    nlat, nlon = size(landcover_mat[1])
    slim, nlim = abs(-60), 80
    res = Int64(nlat/(nlim+slim))
    south_nothings = fill(nothing, ((90-slim)*res, nlon))
    north_nothings = fill(nothing, ((90-nlim)*res, nlon))
    landcover_grids = [vcat(south_nothings, l, north_nothings) for l in landcover_mat]

    # Convert to SimpleSDMLayers
    landcover_layers = SimpleSDMPredictor.(landcover_grids, -180.0, 180.0, -90.0, 90.0)

    # return landcover_layers
# end

[geotiff(SimpleSDMPredictor, p) for p in paths] # always returning 840 x 2159...

coords_landcover = (left = -180.0, right = 180.0, bottom = -60.0, top = 80.0)
geotiff(SimpleSDMPredictor, paths[2]; coords_landcover...) # still the same

## Check details from geotiff
file = paths[2]
bandnumber = 1
left = -180.0
right = 180.0
bottom = -90.0
top = 90.0

# function geotiff(
#     ::Type{LT},
#     file::AbstractString,
#     bandnumber::Integer=1;
#     left = -180.0,
#     right = 180.0,
#     bottom = -90.0,
#     top = 90.0
# ) where {LT<:SimpleSDMLayer}

    # This next block is reading the geotiff file, but also making sure that we
    # clip the file correctly to avoid reading more than we need.
    # This next block is reading the geotiff file, but also making sure that we
    # clip the file correctly to avoid reading more than we need.
    # layer = ArchGDAL.read(file) do dataset
    dataset = ArchGDAL.read(file)

        transform = ArchGDAL.getgeotransform(dataset)
        wkt = ArchGDAL.getproj(dataset)

        # The data we need is pretty much always going to be stored in the first
        # band, but this is not the case for the future WorldClim data.
        band = ArchGDAL.getband(dataset, bandnumber)
        T = ArchGDAL.pixeltype(band)
        
        # The nodata is not always correclty identified, so if it is not found, we assumed it is the smallest value in the band
        nodata = isnothing(ArchGDAL.getnodatavalue(band)) ? convert(T, ArchGDAL.minimum(band)) : convert(T, ArchGDAL.getnodatavalue(band))

        # Get the correct latitudes
        minlon = transform[1]
        maxlat = transform[4]
        maxlon = minlon + size(band,1)*transform[2]
        minlat = maxlat - abs(size(band,2)*transform[6])

        left = isnothing(left) ? minlon : max(left, minlon)
        right = isnothing(right) ? maxlon : min(right, maxlon)
        bottom = isnothing(bottom) ? minlat : max(bottom, minlat)
        top = isnothing(top) ? maxlat : min(top, maxlat)

        lon_stride, lat_stride = transform[2], transform[6]
        
        width = ArchGDAL.width(dataset)
        height = ArchGDAL.height(dataset)

        #global lon_stride, lat_stride
        #global left_pos, right_pos
        #global bottom_pos, top_pos

        lon_stride, left_pos, min_width = SimpleSDMLayers._find_span(width, minlon, maxlon, left)
        _, right_pos, max_width = SimpleSDMLayers._find_span(width, minlon, maxlon, right) # here's the problem!!
        lat_stride, top_pos, max_height = SimpleSDMLayers._find_span(height, minlat, maxlat, top)
        _, bottom_pos, min_height = SimpleSDMLayers._find_span(height, minlat, maxlat, bottom)

        max_height, min_height = height .- (min_height, max_height) .+ 1

        # We are now ready to initialize a matrix of the correct type.
        buffer = Matrix{T}(undef, length(min_width:max_width), length(min_height:max_height))
        ArchGDAL.read!(dataset, buffer, bandnumber, min_height:max_height, min_width:max_width)
        buffer = convert(Matrix{Union{Nothing,eltype(buffer)}}, rotl90(buffer))
        replace!(buffer, nodata => nothing)
        # LT(buffer, left_pos-0.5lon_stride, right_pos+0.5lon_stride, bottom_pos-0.5lat_stride, top_pos+0.5lat_stride)
        SimpleSDMPredictor(buffer, left_pos-0.5lon_stride, right_pos+0.5lon_stride, bottom_pos-0.5lat_stride, top_pos+0.5lat_stride)
    end

    # return layer

# end

geotiff(SimpleSDMPredictor, paths[2]; left=nothing, right=nothing, bottom=nothing, top=nothing)
lc2_geo_full

# Extract problem
_, right_pos, max_width = SimpleSDMLayers._find_span(width, minlon, maxlon, right) # here's the problem!!

# Check details of _find_span
n = width
m = minlon
M = maxlon
pos = right
# function _find_span(n, m, M, pos)
    pos > M && return nothing
    pos < m && return nothing
    stride2 = (M - m) / n
    centers = (m + 0.5stride2):stride2:(M-0.5stride2)
    span_pos = last(findmin(abs.(pos .- centers)))
    return (stride, centers[span_pos], span_pos)
# end

lc2_geo_fix = similar(lc2_csv)
lc2_geo_fix.grid = copy(lc2_geo.grid)
lc2_geo_fix = convert(SimpleSDMPredictor, lc2_geo_fix)

compdf = DataFrame([lc2_geo_fix, lc2_csv])
rename!(compdf, :x1 => :geo_fix, :x2 => :csv)
allowmissing!(compdf)
for col in [:geo_fix, :csv]
    replace!(compdf[!, col], nothing => missing)
end
compdf
dropmissing(compdf, :geo_fix)

plot(lc2_geo_fix, c = :viridis)
plot(lc2_csv, c = :viridis) # one cell misplacement

