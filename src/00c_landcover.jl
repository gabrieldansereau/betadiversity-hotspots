include("required.jl")

## Conditional arguments
# save_prepdata = true # should preparation data be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

## Run bash scripts to download & coarsen landcover data from Zenodo (if files missing)
if (@isdefined save_prepdata) && save_prepdata == true
    # List files in landcover folder
    @info "Starting landcover data saving process"
    lc_files = readdir("assets/landcover/")

    # Check if landcover files are missing
    if !any(startswith.(lc_files, r"^lc_"))
        @info "Landcover files missing. Starting extraction"

        # Check if full resolution files are missing
        if !any(startswith.(lc_files, r"^landcover_copernicus_global_100m"))
            # Download full resolution files
            # BEWARE, can be long to download, 25 GB of data
            @info "Full resolution files missing. Starting download"
            run(`bash src/shell/landcover_download.sh`)
        end

        # Coarsen resolution
        @info "Coarsening resolution"
        run(`bash src/shell/landcover_coarsen.sh`)
    end
else
    @info "Not saving landcover data"
end

## Test landcover variables
# Define coordinates range
coords = (left=-145.0, right=-50.0, bottom=20.0, top=75.0)

# Test loading variables
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution=5.0)
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution=10.0)
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution=10.0, coords...)

# Load worldclim variables to compare
wc_vars = SimpleSDMPredictor(WorldClim, BioClim, 1:19; resolution=10.0, coords...);

## Export to CSV

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

# Create env dataframe (environmental data)
env_df = DataFrame(convert.(Float64, env_vars))
select!(env_df, Not([:longitude, :latitude]))
rename!(env_df, vcat(Symbol.("wc", 1:size(wc_vars, 1)), Symbol.("lc", 1:size(lc_vars, 1))))
insertcols!(env_df, 1, :site => 1:nrow(env_df))

# Create spa DataFrame (spatial coordinates)
spa_df = DataFrame(wc_vars[1])
select!(spa_df, Not(:values))
insertcols!(spa_df, 1, :site => Float64.(1:nrow(spa_df)))
rename!(spa_df, [:site, :lon, :lat])

# Create spa layers
spa_vars = [
    SimpleSDMPredictor(spa_df, x, wc_vars[1]; latitude=:lat, longitude=:lon) for
    x in [:site, :lon, :lat]
]

# Export dataframes
# save_prepdata = true
if (@isdefined save_prepdata) && save_prepdata == true
    @info "Exporting env & spa to CSV"
    CSV.write(joinpath("data", "proc", "distributions_env_full.csv"), env_df; delim="\t")
    CSV.write(joinpath("data", "proc", "distributions_spa_full.csv"), spa_df; delim="\t")
    geotiff(joinpath("data", "raster", "env_stack.tif"), env_vars)
    geotiff(joinpath("data", "raster", "spa_stack.tif"), spa_vars)
end

# Test load
#=
testspa = CSV.read(joinpath("data", "proc", "distributions_spa_full.csv"), DataFrame, header=true, delim="\t")
testenv = CSV.read(joinpath("data", "proc", "distributions_env_full.csv"), DataFrame, header=true, delim="\t")
testjoin = innerjoin(testspa, testenv, on = :site) =#

## Export QC sites coordinates (for smaller scale analyses)

# Define Quebec coordinates
coords_qc = (left=-80.0, right=-55.0, bottom=45.0, top=63.0)

# Get site indices
spa_qc = filter(:lon => >=(coords_qc.left - stride(wc_vars[1]; dims=2)), spa_df)
filter!(:lon => <(coords_qc.right), spa_qc)
filter!(:lat => >=(coords_qc.bottom - stride(wc_vars[1], dims=1)), spa_qc)
filter!(:lat => <(coords_qc.top), spa_qc)

# Subset layers to QC only
env_vars_qc = [v[coords_qc] for v in env_vars]
spa_vars_qc = [v[coords_qc] for v in spa_vars]

# Export QC dataframes
# save_prepdata = true
if (@isdefined save_prepdata) && save_prepdata == true
    @info "Exporting QC env & spa to CSV"
    CSV.write(joinpath("data", "proc", "distributions_spa_qc.csv"), spa_qc; delim="\t")
    geotiff(joinpath("data", "raster", "env_stack_qc.tif"), env_vars_qc)
    geotiff(joinpath("data", "raster", "spa_stack_qc.tif"), spa_vars_qc)
end
