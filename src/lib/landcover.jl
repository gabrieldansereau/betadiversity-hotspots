## Landcover variables loading functions

# Landcover variables loading
function landcover(layers::Vector{Int64}; resolution::AbstractString="10", path::AbstractString="assets/landcover/")
    ## Get file paths
    # List files in path
    lc_files = readdir("$(path)")
    # Filter for selected resolution
    filter!(x -> occursin.("$(resolution)m.tif", x), lc_files)
    # Create paths for selected layers only
    paths = joinpath.(path, lc_files)[layers]

    ## Load data
    # Read raster dataset
    datasets = [ArchGDAL.read(p) for p in paths]
    # Read raster values from band 1 (only band in this case)
    values = [ArchGDAL.read(d, 1) for d in datasets]
    # Convert from UInt8 to Float64
    values_int = [convert.(Float64, v) for v in values]

    ## Preliminary manipulations
    # Reverse rows and permute dimensions
    landcover_mat = [permutedims(v[:,end:-1:1]) for v in values_int]
    # Replace 255 (default no data values) by NaN
    [replace!(l, 255 => NaN) for l in landcover_mat]
    # Add additionnal line of NaNs down South (coordinates prime did not exactly match)
    landcover_mat = [vcat(fill(NaN, (1, size(l, 2))), l) for l in landcover_mat]
    # Convert to SDMLayers
    landcover_layers = SimpleSDMPredictor.(landcover_mat, -160.0, -40.0, 20.0, 80.0)
end

"""
    landcover(layer::Int64; x...)

Return a single landcover layer.
"""
landcover(layer::Int64; x...) = landcover([layer]; x...)[1]

"""
    landcover(layers::UnitRange{Int64}; x...)

Return a range of landcover layers.
"""
landcover(layers::UnitRange{Int64}; x...) = landcover(collect(layers); x...)
