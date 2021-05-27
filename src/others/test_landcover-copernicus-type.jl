# struct Copernicus end

# function SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layer::Integer=1; kwargs...)
#     file = _get_raster(EarthEnv, LandCover, layer, full)
#     return geotiff(SimpleSDMPredictor, file; kwargs...)
# end

# function SimpleSDMPredictor(::Type{EarthEnv}, ::Type{LandCover}, layers::AbstractArray; kwargs...)
#     @assert eltype(layers) <: Integer
#     return [SimpleSDMPredictor(EarthEnv, LandCover, l; kwargs...) for l in layers]
# end

import Pkg; Pkg.activate(".")
include("../required.jl")

coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)

# Test loading variables
lc_vars = landcover(1:10, resolution = 5.0)
lc_vars = landcover(1:10, resolution = 10.0)
lc_vars = map(x -> landcover(x, resolution = 10.0)[coords], 1:10)

# Test loading variables
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution = 5.0)
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10)
lc_vars = map(x -> SimpleSDMPredictor(Copernicus, LandCover, x)[coords], 1:10)

# Check compatibility
lc1 = SimpleSDMPredictor(Copernicus, LandCover, 1)
wc1 = SimpleSDMPredictor(WorldClim, BioClim, 1)
SimpleSDMLayers._layers_are_compatible(lc1, wc1) # compatible

# lc1 = SimpleSDMPredictor(Copernicus, LandCover, 1; coords...)
lc1 = SimpleSDMPredictor(Copernicus, LandCover, 1)[coords]
lc1 = lc1[top = coords.top - stride(lc1, 2)]
wc1 = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)
SimpleSDMLayers._layers_are_compatible(lc1, wc1) # not compatible