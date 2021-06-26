include("required.jl")

## Conditional arguments
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
# save_data = true # should data files be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

## Get & prepare data

# Define coordinates range
coords = (left=-145.0, right=-50.0, bottom=20.0, top=75.0)
copy_layer = SimpleSDMPredictor(WorldClim, BioClim, 1; coords...)

# Load data from CSV files
df = CSV.read(joinpath("data", "proc", "ebd_warblers_prep.csv"), DataFrame)

# Remove groupIdentifier column (causing bug)
select!(df, Not(:groupIdentifier))

# Filter observations outside coordinates range
filter!(:longitude => >(coords.left), df)
filter!(:longitude => <(coords.right), df)
filter!(:latitude => >(coords.bottom), df)
filter!(:latitude => <(coords.top), df)

# Separate warblers species
warblers = groupby(df, :species)

# Reorder species by frequency
warblers = warblers |>
    x -> combine(nrow, x) |>
    x -> sortperm(x, :nrow, rev=true) |>
    x -> warblers[x]

# Extract species names
spenames = [w.species[1] for w in warblers]
specommon = [w.commonName[1] for w in warblers]

# Create index Dict for species names
speindex = indexmap(spenames)

## Get environmental data

# WorldClim data
wc_vars = SimpleSDMPredictor(WorldClim, BioClim, [1, 12]; resolution=10.0, coords...)

# Landcover data
lc_vars = SimpleSDMPredictor(Copernicus, LandCover, 1:10; resolution=10.0, coords...)

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)

## Get distribution for all species
# Runs only if create_distributions = true, as it can take a while. Loaded from file otherwise
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
if (@isdefined create_distributions) && create_distributions == true
    @info "Raw distributions to be created"
    # Get raw distributions
    @time distributions = @showprogress [presence_absence(w, env_vars[1]) for w in warblers]
end

## Export distributions
# save_data = true # should data files be overwritten (optional)
if (@isdefined save_data) && save_data == true
    ## Export to JLD2
    # Export data to JLD2
    @info "Exporting data to JLD2 file (raw distributions data)"
    @save joinpath("data", "jld2", "raw-distributions.jld2") distributions spenames specommon speindex
    # Export to ZIP archive
    @info "Exporting data from JLD2 to ZIP archive (raw distributions data)"
    _zip_jld2(
        joinpath("data", "jld2", "raw-distributions.zip"),
        joinpath("data", "jld2", "raw-distributions.jld2"),
    )
    # Make sure JLD2 timestamp is more recent than ZIP archive
    touch(joinpath("data", "jld2", "raw-distributions.jld2"))

    ## Export species names
    @info "Exporting species names to JLD2 file"
    @save joinpath("data", "jld2", "spenames.jld2") spenames specommon speindex

    ## Export glossary
    @info "Exporting species names to glossary"
    # Add species names in correct order to glossary
    include(joinpath("others", "data_glossary.jl"))

    ## Export layers
    # Export distribution layers
    geotiff(joinpath("data", "proc", "distributions_raw.tif"), distributions)
    # Export distributions for QC only
    coords_qc = (left=-80.0, right=-55.0, bottom=45.0, top=63.0)
    distributions_qc = [d[coords_qc] for d in distributions]
    geotiff(joinpath("data", "proc", "distributions_raw_qc.tif"), distributions_qc)

    ## Export to CSV as Y matrix
    @info "Exporting data to CSV as Y matrix (raw distributions data)"
    # Get Y matrix
    Y = calculate_Y(distributions)
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)
    # Convert to dataframe
    spe_df = DataFrame(Yobs, :auto)
    rename!(spe_df, Symbol.("sp", 1:ncol(spe_df)))
    insertcols!(spe_df, 1, :site => inds_obs)
    # Export to CSV
    CSV.write(joinpath("data", "proc", "distributions_spe_full.csv"), spe_df; delim="\t")
else
    # Load data
    @info "Data imported from file (raw distributions data)"
    @load joinpath("data", "jld2", "raw-distributions.jld2") spenames specommon speindex
    distributions = [
        geotiff(
            SimpleSDMPredictor, joinpath("data", "proc", "distributions_raw.tif"), i
        ) for i in eachindex(spenames)
    ]
end

## Count sites with presence per species
pres_counts = [
    length(filter(x -> !isnothing(x) && x > 0.0, species.grid)) for species in distributions
]
sort(pres_counts)

## Plot result
# Species 1
sp1 = "Setophaga_townsendi"
map_sp1 = plotSDM2(
    distributions[speindex[sp1]];
    c=:BuPu,
    title="$(replace(sp1, "_" => " ")) distribution (raw)",
    colorbar=:none,
    dpi=200,
)
# Species 2
sp2 = "Setophaga_petechia"
map_sp2 = plotSDM2(
    distributions[speindex[sp2]];
    c=:BuPu,
    title="$(replace(sp2, "_" => " ")) distribution (raw)",
    colorbar=:none,
    dpi=200,
)
