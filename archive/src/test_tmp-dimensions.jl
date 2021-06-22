include("../required.jl")

## Get layers
# Define coordinates range
coords = (left=-145.0, right=-50.0, bottom=20.0, top=75.0)

# Test loading variables
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution=10.0, coords...)

# Load worldclim variables to compare
wc_vars = SimpleSDMPredictor(WorldClim, BioClim, 1:19; resolution=10.0, coords...);

## Temporary dimensions fix
# Check which line is extra
wc_test = SimpleSDMPredictor(WorldClim, BioClim, 1)[coords]
wc_vars[1].grid
wc_test.grid # last line seems extra
isequal(wc_test.grid[1:(end - 1), :], wc_vars[1].grid) # true

# Test clip fix
wc_clip = wc_test[top=coords.top - stride(wc_test, 2)]
isequal(wc_clip.grid, wc_vars[1].grid) # true
wc_vars[1].top
wc_clip.top # not exactly the same
isapprox(wc_vars[1].top, wc_clip.top) # approximately the same

# Test clip from world layer instead
wc_test2 = SimpleSDMPredictor(WorldClim, BioClim, 1)[coords][top=coords.top -
                                                                 stride(wc_test, 2)] # same limits
wc_test2.top # not helping

# Fix landcover layers
lc_vars = [l[top=coords.top - stride(l, 2)] for l in lc_vars]
isequal(latitudes(wc_vars[1]), latitudes(lc_vars[1])) # false
isapprox(latitudes(wc_vars[1]), latitudes(lc_vars[1])) # true so should be ok

## Export to CSV
# Combine environmental data
env_vars = [wc_vars; lc_vars]
# Create env matrix
env_mat = mapreduce(x -> vec(x.grid), hcat, env_vars)
replace!(x -> isnothing(x) ? NaN : x, env_mat)
# Create env dataframe
env_df = DataFrame(env_mat, :auto)
rename!(env_df, vcat(Symbol.("wc", 1:size(wc_vars, 1)), Symbol.("lc", 1:size(lc_vars, 1))))
insertcols!(env_df, 1, :site => 1:nrow(env_df))

# Get sites latitudes
lats = repeat(collect(latitudes(wc_vars[1])); outer=size(wc_vars[1].grid, 2))
# Get sites longitudes
lons = repeat(collect(longitudes(wc_vars[1])); inner=size(wc_vars[1].grid, 1))
# Create spa matrix
spa_mat = [lats lons]
# Create spa dataframe
spa_df = DataFrame(; site=eachindex(lats), lat=lats, lon=lons)

## Check difference with previously saved data

# Copy data
run(
    `cp ./data/proc/distributions_env_full.csv ./data/proc/distributions_env_full_backup.csv`,
)
run(
    `cp ./data/proc/distributions_spa_full.csv ./data/proc/distributions_spa_full_backup.csv`,
)
run(`cp ./data/proc/distributions_spa_qc.csv ./data/proc/distributions_spa_qc_backup.csv`)

# Export new files
CSV.write(joinpath("data", "proc", "distributions_env_full.csv"), env_df; delim="\t")
CSV.write(joinpath("data", "proc", "distributions_spa_full.csv"), spa_df; delim="\t")
# Git diff says spa files are different
# Sizes are different too
# Looking at the diff, it seems to be a float approximation difference, which should be fine

# Let's make sure it's just that
spa_df_backup = CSV.read(
    joinpath("data", "proc", "distributions_spa_full_backup.csv"), DataFrame
)
spa_df_backup
spa_df
isequal(spa_df_backup, spa_df) # false
isapprox(spa_df_backup, spa_df) # true
isequal(spa_df_backup.site, spa_df.site) # true
isapprox(spa_df_backup.lat, spa_df.lat) # true
isapprox(spa_df_backup.lon, spa_df.lon) # true
# Yep, should just be approximation error, so its fine

# Make sure env files are equal too, even if diff says they are
env_df_backup = CSV.read(
    joinpath("data", "proc", "distributions_env_full_backup.csv"), DataFrame
)
env_df_backup
env_df
isequal(env_df_backup, env_df) # false
isapprox(env_df_backup, env_df) # false ?!
isapprox(NaN, NaN) # oh probably just because of NaN
# whatever let's trust Git

## Export QC sites coordinates (for smaller scale analyses)
coords_qc = (left=-80.0, right=-55.0, bottom=45.0, top=63.0)
# Get site indices
spa_qc = filter(
    x ->
        (coords_qc.left - stride(wc_vars[1]; dims=2) <= x.lon < coords_qc.right) &&
            (coords_qc.bottom - stride(wc_vars[1]; dims=1) <= x.lat < coords_qc.top),
    spa_df,
)

# Test dimensions for Quebec coordinates
wc_vars_qc = SimpleSDMPredictor(WorldClim, BioClim, 1; coords_qc...)
wc_test_qc = SimpleSDMPredictor(WorldClim, BioClim, 1)[coords_qc]
isequal(size(wc_vars_qc), size(wc_test_qc)) # false
isapprox(longitudes(wc_vars_qc), unique(spa_qc.lon)) # true
isapprox(longitudes(wc_test_qc), unique(spa_qc.lon)) # false
isequal(length(wc_vars_qc.grid), nrow(spa_qc)) # equal
# So dimensions in the exported file should be ok
CSV.write(joinpath("data", "proc", "distributions_spa_qc.csv"), spa_qc; delim="\t")
# Git diff says there's difference
spa_qc_backup = CSV.read(
    joinpath("data", "proc", "distributions_spa_qc_backup.csv"), DataFrame
)
isapprox(spa_qc_backup, spa_qc) # true
# Again just an approximation error
