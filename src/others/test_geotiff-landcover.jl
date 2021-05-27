import Pkg; Pkg.activate(".")
include("../required.jl")

## Test landcover variables
# Define coordinates range
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)

# Test loading variables
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution = 10.0, coords...)

# Load worldclim variables to compare
wc_vars = SimpleSDMPredictor(WorldClim, BioClim, 1:19; resolution = 10.0, coords...)

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

## Attempt to load full 100m resolution layer
megafilepath = "./assets/landcover/test/landcover_copernicus_global_100m_v2.0.2_moss.tif"
megafilepath = "./assets/landcover/landcover_copernicus_global_100m_v2.0.2_urban.tif"
isfile(megafilepath)
test = geotiff(SimpleSDMPredictor, megafilepath; coords...)
ArchGDAL.read(megafilepath)
dataset = ArchGDAL.read(megafilepath)
ArchGDAL.getgeotransform(dataset)
ArchGDAL.read("./assets/landcover/lc_bare_10m.tif")
ArchGDAL.getgeotransform(ArchGDAL.read("./assets/landcover/lc_bare_10m.tif"))

wc1_path = "/home/gdansereau/.data/SimpleSDMLayers.jl/assets/WorldClim/BioClim/10/wc2.1_10m_bio_1.tif"
geotiff(SimpleSDMPredictor, wc1_path)
wc1_dataset = ArchGDAL.read(wc1_path)
ArchGDAL.getgeotransform(wc1_dataset)

## In bash
# -tap option
gdalwarp -tr 0.166667 0.166667 -tap -r average --config GDAL_CACHEMAX 500 -wm 500 -multi assets/landcover/lc_moss_10m.tif -overwrite assets/landcover/landcover_test.tif

# -te option
gdalwarp -tr 0.166667 0.166667 -te -180.0 -60.0 180.0 80.0 -r average --config GDAL_CACHEMAX 500 -wm 500 -multi assets/landcover/lc_moss_10m.tif -overwrite assets/landcover/landcover_test.tif

# -te option and same stride as WorldClim
gdalwarp -tr 0.1666666666666666574 0.1666666666666666574 -te -180.0 -60.0 180.0 80.0 -r average --config GDAL_CACHEMAX 500 -wm 500 -multi assets/landcover/lc_moss_10m.tif -overwrite assets/landcover/landcover_test.tif

# -te option and same stride as WorldClim on full layer
gdalwarp -tr 0.1666666666666666574 0.1666666666666666574 -te -180.0 -60.0 180.0 80.0 -r average --config GDAL_CACHEMAX 500 -wm 500 -multi assets/landcover/landcover_copernicus_global_100m_v2.0.2_moss.tif -overwrite assets/landcover/landcover_fulltest.tif

# full spatial extent
gdalwarp -tr 0.1666666666666666574 0.1666666666666666574 -te -180.0 -90.0 180.0 90.0 -r average --config GDAL_CACHEMAX 500 -wm 500 -multi assets/landcover/landcover_copernicus_global_100m_v2.0.2_moss.tif -overwrite assets/landcover/landcover_fulltest.tif

# Back in Julia
test_path = "./assets/landcover/landcover_test.tif"
test_dataset = ArchGDAL.read(test_path)
ArchGDAL.getgeotransform(test_dataset)

moss_path = "./assets/landcover/lc_moss_10m.tif"
moss_dataset = ArchGDAL.read(moss_path)
ArchGDAL.getgeotransform(moss_dataset)

fulltest_path = "./assets/landcover/landcover_fulltest.tif"
fulltest_dataset = ArchGDAL.read(fulltest_path)
ArchGDAL.getgeotransform(fulltest_dataset)

geotiff(SimpleSDMPredictor, test_path)
geotiff(SimpleSDMPredictor, fulltest_path)
test_layer = convert(Float32, geotiff(SimpleSDMPredictor, test_path; coords...))
fulltest_layer = convert(Float32, geotiff(SimpleSDMPredictor, fulltest_path; coords...))
wctest_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
moss_indx = readdir("assets/landcover/") |> x -> 
    filter(f -> occursin.("10m.tif", f), x) |> x ->
    findfirst(occursin.("moss", x))
lc_moss = lc_vars[moss_indx]

isequal(latitudes(fulltest_layer), latitudes(wctest_layer))
isequal(longitudes(fulltest_layer), longitudes(wctest_layer))
SimpleSDMLayers._layers_are_compatible(fulltest_layer, wctest_layer)

plot(test_layer)
plot(lc_moss)
plot(fulltest_layer)

isequal(lc_moss, fulltest_layer)
isequal(lc_moss.grid, fulltest_layer.grid)
testdf = DataFrame([lc_moss, fulltest_layer]) # not exactly same coordinates

isapprox(lc_moss.bottom, fulltest_layer.bottom)

moss_new = similar(lc_moss)
moss_new.grid = copy(fulltest_layer.grid)
moss_new = convert(SimpleSDMPredictor, moss_new)

testdf = DataFrame([lc_moss, moss_new])
isequal(testdf.x1, testdf.x2)
isapprox(testdf.x1, testdf.x2)
isequal(collect(lc_moss), collect(moss_new))
isapprox(collect(lc_moss), collect(moss_new))

writedlm("test.csv", collect(lc_moss))
writedlm("test.csv", collect(moss_new))

allowmissing!(testdf)
for col in [:x1, :x2]
    replace!(testdf[!, col], nothing => missing)
end
testdf
filter(x -> ismissing(x.x1) && !ismissing(x.x2), testdf)
filter(x -> !ismissing(x.x1) && ismissing(x.x2), testdf)
isequal(nrow(dropmissing(testdf)), length(lc_moss))

dropmissing!(testdf)
testdiff = filter(x -> x.x1 != x.x2, testdf) # 800 sites with difference
insertcols!(testdiff, :diff => testdiff.x2 .- testdiff.x1)
unique(testdiff.diff)
sort(countmap(testdiff.diff))