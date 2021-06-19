## Landcover variables loading functions

struct Copernicus end

import SimpleSDMLayers: SimpleSDMPredictor

"""
    SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layer::Integer=1; resolution::Float64=10.0, path::AbstractString=joinpath("assets", "landcover"), kwargs...)
    SimpleSDMPredictor(::Type{Copernicus}, ::Type{LandCover}, layers::AbstractArray; kwargs...)

Prepares the Copernicus v2.0.1 landcover data, and returns them as an array of 
`SimpleSDMPredictor`s. Layers are called by their number, from 1 to 10. The list
of available layers is given in the table below.

This function loads the pre-downloaded and coarsened data from the 
`assets/landcover/` folder (the `path` argument). These steps can be performed 
with the scripts in `src/shell/` using Bash and GDAL. The available resolutions 
are 10 and 5 arc-minutes, and must be specified as a floating point value with 
the `resolution` keyword (with `10.0` as the default).

The original data by Buchhorn et al. is archived on Zenodo at 
https://zenodo.org/record/3243509, and is available at a much finer resolution.

| Variable | Description        |
| ------   | ------             |
| 1        | Bare               |
| 2        | Crops              |
| 3        | Grass              |
| 4        | Moss               |
| 5        | Shrub              |
| 6        | Snow               |
| 7        | Tree               |
| 8        | Urban              |
| 9        | Water permanent    |
| 10       | Water seasonal     |

"""
function SimpleSDMPredictor(
    ::Type{Copernicus},
    ::Type{LandCover},
    layer::Integer=1;
    resolution::Float64=10.0,
    path::AbstractString=joinpath("assets", "landcover"),
    kwargs...,
)
    @assert resolution in [5.0, 10.0]
    1 ≤ layer ≤ 10 || throw(ArgumentError("The layer must be between 1 and 10"))

    # List files in path
    lc_files = readdir(path)
    # Filter for selected resolution
    filter!(contains("$(Int(resolution))m.tif"), lc_files)
    # Create path for selected layer only
    p = joinpath.(path, lc_files)[layer]

    # Read layer
    layer = geotiff(SimpleSDMPredictor, p; kwargs...)
    # Convert from UInt8 to Float32
    layer = convert(Float32, layer)

    return convert(SimpleSDMPredictor, layer)
end

function SimpleSDMPredictor(
    ::Type{Copernicus}, ::Type{LandCover}, layers::AbstractArray; kwargs...
)
    @assert eltype(layers) <: Integer
    return [SimpleSDMPredictor(Copernicus, LandCover, l; kwargs...) for l in layers]
end
