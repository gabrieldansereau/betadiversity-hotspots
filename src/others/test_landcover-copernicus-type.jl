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

lc_vars = map(x -> SimpleSDMPredictor(Copernicus, LandCover, x)[coords], 1:10)
lc_vars = [l[top = coords.top - stride(l, 2)] for l in lc_vars]

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
