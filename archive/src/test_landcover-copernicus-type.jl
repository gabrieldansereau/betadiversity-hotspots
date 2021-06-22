import Pkg; Pkg.activate(".")
include("../required.jl")

coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)

# Test loading variables
# lc_vars = landcover(1:10, resolution = 5.0)
# lc_vars = landcover(1:10, resolution = 10.0)
# lc_vars = map(x -> landcover(x, resolution = 10.0)[coords], 1:10)

# # Test loading variables
# lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution = 5.0)
# lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10)
# lc_vars = map(x -> SimpleSDMPredictor(Copernicus, LandCover, x)[coords], 1:10)

# # Check compatibility
# lc1 = SimpleSDMPredictor(Copernicus, LandCover, 1)
# wc1 = SimpleSDMPredictor(WorldClim, BioClim, 1)
# SimpleSDMLayers._layers_are_compatible(lc1, wc1) # compatible

# # lc1 = SimpleSDMPredictor(Copernicus, LandCover, 1; coords...)
# lc1 = SimpleSDMPredictor(Copernicus, LandCover, 1)[coords]
# lc1 = lc1[top = coords.top - stride(lc1, 2)]
# wc1 = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
# SimpleSDMLayers._layers_are_compatible(lc1, wc1) # not compatible

## Export values to compare
lc_path = joinpath("assets", "landcover")
function nothingtomissing!(df::DataFrame)
    allowmissing!(df)
    for col in eachcol(df)
        replace!(col, nothing => missing)
    end
end
function valuesdf(layers::Vector{T}) where {T <: SimpleSDMLayer}
    df = DataFrame(layers)
    select!(df, Not([:longitude, :latitude]))
    rename!(df, Symbol.("lc", 1:10))
    nothingtomissing!(df)
    return df
end

# lc_vars = map(x -> SimpleSDMPredictor(Copernicus, LandCover, x)[coords], 1:10)
# lc_vars = [l[top = coords.top - stride(l, 2)] for l in lc_vars]
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10)
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; coords...)
lc1 = copy(lc_vars)[1]
wc1 = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)

# Step 1: SimpleSDMPredictor, as before
df1 = valuesdf(lc_vars)
# CSV.write(joinpath(lc_path, "landcover_values_df1.csv"), df1)
df1 = CSV.read(joinpath(lc_path, "landcover_values_df1.csv"), DataFrame)

# Step 2: SimpleSDMPredictor, with integer as default and range as overload
SimpleSDMPredictor(Copernicus, LandCover)
SimpleSDMPredictor(Copernicus, LandCover, 1)
SimpleSDMPredictor(Copernicus, LandCover, 1:10)
df2 = valuesdf(lc_vars)
# CSV.write(joinpath(lc_path, "landcover_values_df2.csv"), df2)
df2 = CSV.read(joinpath(lc_path, "landcover_values_df2.csv"), DataFrame)
isequal(df1, df2) # true
# confirmed by comparison

# Step 3. load with geotiff instead of custom ArchGDAL call
layer = 1
resolution = 10.0
path = joinpath("assets", "landcover")

# List files in path
lc_files = readdir(path)
# Filter for selected resolution
filter!(x -> occursin.("$(Int(resolution))m.tif", x), lc_files)
# Create path for selected layer only
p = joinpath.(path, lc_files)[layer]
# Use geotiff instead
geotiff(SimpleSDMPredictor, p)
ArchGDAL.getgeotransform(d)
wctransform = ArchGDAL.getgeotransform(ArchGDAL.read(joinpath(ENV["SDMLAYERS_PATH"], "WorldClim", "BioClim", "10", "wc2.1_10m_bio_1.tif")))

gdalwarp -te -180.0 -60.0 180.0 80.0 assets/landcover/lc_bare_10m.tif assets/landcover/lc_bare_10m2.tif

p2 = replace(p, ".tif" => "2.tif")
d2 = ArchGDAL.read(p2)
ArchGDAL.getgeotransform(d2)
landcover_layers2 = geotiff(SimpleSDMPredictor, p2)

l3 = geotiff(SimpleSDMPredictor, p2; coords...)
SimpleSDMLayers._layers_are_compatible(l3, l1) # not compatible
SimpleSDMLayers._layers_are_compatible(l3, wc1) # true, compatible with WorldClim!!
isequal(l3.grid, l1.grid) # grids are equal!!

# Export values
lc_files = readdir(path)
filter!(x -> occursin.("$(Int(resolution))m.tif", x), lc_files)
ps = joinpath.(path, lc_files)
lc_vars = geotiff.(SimpleSDMPredictor, ps)
lc_vars = geotiff.(SimpleSDMPredictor, ps; coords...)
lc_vars = convert.(Float32, lc_vars)
df3 = valuesdf(lc_vars)
# CSV.write(joinpath(lc_path, "landcover_values_df1.csv"), df3)
# CSV.write(joinpath(lc_path, "landcover_values_df3.csv"), df3)
df3 = CSV.read(joinpath(lc_path, "landcover_values_df3.csv"), DataFrame)
isequal(df1, df3) # true
# confirmed by comparison

# Step 4: Update function
SimpleSDMLayers._layers_are_compatible(lc1, wc1)
df4 = valuesdf(lc_vars)
# CSV.write(joinpath(lc_path, "landcover_values_df1.csv"), df4)
# CSV.write(joinpath(lc_path, "landcover_values_df4.csv"), df4)
df4 = CSV.read(joinpath(lc_path, "landcover_values_df4.csv"), DataFrame)
isequal(df1, df4) # true
# confirmed by comparison

# Step 5: stack
stack_path = joinpath("assets", "landcover", "landcover.tif")
geotiff(stack_path, lc_vars)
lc_vars2 = [geotiff(SimpleSDMPredictor, stack_path, i) for i in 1:10]
lc_vars2 = replace.(lc_vars2, NaN => nothing)
isequal(lc_vars, lc_vars2)
isequal(lc_vars[1], lc_vars2[1])
isequal(lc_vars[1].grid, lc_vars2[1].grid)
isequal(size(lc_vars[1]), size(lc_vars2[1].grid)) # not same size...
# so dimensions problem with geotiff writing and/or reading call too...

geotiff("test.tif", lc_vars[1])
test1 = geotiff(SimpleSDMPredictor, "test.tif")
isequal(size(lc_vars[1]), test1) # not equal
dtest1 = ArchGDAL.read("test.tif") # dimensions are good in tif file
ArchGDAL.getgeotransform(dtest1) # possibly are problem with left coordinate?