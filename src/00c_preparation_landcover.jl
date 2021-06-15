import Pkg; Pkg.activate(".")
include("required.jl")

## Conditional arguments
# save_prepdata = true # should preparation data be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

## Run bash scripts to download & coarsen landcover data from Zenodo (if files missing)
if (@isdefined save_prepdata) && save_prepdata == true
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
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)

# Test loading variables
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution = 5.0)
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution = 10.0)
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution = 10.0, coords...)

# Load worldclim variables to compare
wc_vars = SimpleSDMPredictor(WorldClim, BioClim, 1:19; resolution = 10.0, coords...);

## Plot environmental variables examples
# Plot wcvars1 (temperature)
wc_plot = plotSDM2(wc_vars[1], colorbar_title = "Annual Mean Temperature (°C)")
# Plot lcvars2 (urban)
lc_plot = plotSDM2(lc_vars[2], colorbar_title = "Crops land cover (%)")

# Get variable info
glossary = CSV.read(joinpath("data", "proc", "glossary.csv"), DataFrame)
lcdf, wcdf, spdf = groupby(glossary, :type)

# Clim variables
heatplots = [plotSDM2(wc_vars[i], 
                      c = :OrRd, 
                      title = uppercasefirst(replace(wcdf.full_name[i], "_" => " ")), 
                      colorbar_title = wcdf.description[i])
            for i in (1,2,5,6)]
heatplot = plot(heatplots..., size = (900, 900))
precplots = [plotSDM2(wc_vars[i], 
                      c = :PuBu, 
                      title = uppercasefirst(replace(wcdf.full_name[i], "_" => " ")), 
                      colorbar_title = wcdf.description[i])
            for i in (12,13,14,15)]
precplot = plot(precplots..., size = (900, 900))

# Landcover variables
landplots = [plotSDM2(lc_vars[i],
                     c = :BuGn,
                     title = uppercasefirst(replace(lcdf.full_name[i], "_" => " ")), 
                     colorbar_title = lcdf.description[i])
             for i in (1:5..., 7:10...)]
landplot = plot(landplots..., size = (1200, 900))


## Export figures
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved (environmental variables)"
    savefig(plot(wc_plot, dpi = 300), joinpath("fig", "00c_wc1-temperature.png"))
    savefig(plot(lc_plot, dpi = 300), joinpath("fig", "00c_lc2-crops.png"))
    savefig(plot(heatplot, dpi = 300), joinpath("fig", "00c_wc-temperatures.png"))
    savefig(plot(precplot, dpi = 300), joinpath("fig", "00c_wc-precipitations.png"))
    savefig(plot(landplot, dpi = 300), joinpath("fig", "00c_lc-landcovers.png"))
else
    @info "Figures not saved (environmental variables)"
end

## Export to CSV

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)
# Create env matrix
env_mat = mapreduce(x -> vec(x.grid), hcat, env_vars)
replace!(env_mat, nothing => NaN)
# Create env dataframe
env_df = DataFrame(env_mat, :auto)
rename!(env_df, vcat(Symbol.("wc", 1:size(wc_vars, 1)), Symbol.("lc", 1:size(lc_vars, 1))))
insertcols!(env_df, 1, :site => 1:nrow(env_df))

# Create spa DataFrame
spa_df = DataFrame(wc_vars[1])
select!(spa_df, Not(:values))
insertcols!(spa_df, 1, :site => Float64.(1:nrow(spa_df)))
rename!(spa_df, [:site, :lon, :lat])
# Create spa layers
spa_vars = [SimpleSDMPredictor(spa_df, x, wc_vars[1]; latitude=:lat, longitude=:lon) for x in [:site, :lon, :lat]]

# Export dataframes
# save_prepdata = true
if (@isdefined save_prepdata) && save_prepdata == true
    @info "Exporting env & spa to CSV"
    CSV.write(joinpath("data", "proc", "distributions_env_full.csv"), env_df, delim="\t")
    CSV.write(joinpath("data", "proc", "distributions_spa_full.csv"), spa_df, delim="\t")
    geotiff(joinpath("data", "proc", "env_stack.tiff"), env_vars)
    geotiff(joinpath("data", "proc", "spa_stack.tiff"), spa_vars)
end

# Test load
#=
testspa = CSV.read(joinpath("data", "proc", "distributions_spa_full.csv"), DataFrame, header=true, delim="\t")
testenv = CSV.read(joinpath("data", "proc", "distributions_env_full.csv"), DataFrame, header=true, delim="\t")
testjoin = innerjoin(testspa, testenv, on = :site)
=#

## Export QC sites coordinates (for smaller scale analyses)

coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
# Get site indices
spa_qc = filter(x -> (coords_qc.left - stride(wc_vars[1], dims = 2) <= x.lon < coords_qc.right) &&
                     (coords_qc.bottom - stride(wc_vars[1], dims = 1) <= x.lat < coords_qc.top), spa_df)

# Export QC dataframes
# save_prepdata = true
if (@isdefined save_prepdata) && save_prepdata == true
    @info "Exporting QC env & spa to CSV"
    CSV.write(joinpath("data", "proc", "distributions_spa_qc.csv"), spa_qc, delim="\t")
end
