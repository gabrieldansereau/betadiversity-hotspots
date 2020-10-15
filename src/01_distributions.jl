if !(@isdefined BetadiversityHotspots)
    import Pkg; Pkg.activate(".")
    using Distributed
    addprocs(9)
    @time include("required.jl")
end

## Conditional arguments
outcome = "raw" # desired outcome (required)
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
# save_data = true # should data files be overwritten (optional)
# save_figures = true # should figures be overwritten (optional)

# Make sure "outcome" is defined
if !(@isdefined outcome)
    @warn "'outcome' not defined, must be either 'raw' or 'bio'"
elseif (outcome != "raw" && outcome != "bio")
    @warn "'outcome' invalid, must be either 'raw' or 'bio'"
else
    @info "'outcome' currently set to '$(outcome)'"
end

## Get & prepare data
# Define coordinates range
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
# Load data from CSV files
df = DataFrame!(CSV.File(joinpath("data", "proc", "ebd_warblers_prep.csv"), header=true, delim="\t"))
# Remove groupIdentifier column (causing bug)
select!(df, Not(:groupIdentifier))
# Filter observations outside coordinates range
filter!(x -> coords.left < x.longitude < coords.right, df)
filter!(x -> coords.bottom < x.latitude < coords.top, df)
# Separate species
# warblers = [df[df.species .== u,:] for u in unique(df.species)]
warblers = groupby(df, :species)
# Reorder species by frequency
# sort!(warblers, by = x -> nrow(x), rev = true)
warblers = warblers |> 
    x -> combine(nrow, x) |>
    x -> sortperm(x, :nrow, rev = true) |>
    x -> warblers[x]
# Extract species names
spenames = [w.species[1] for w in warblers]
specommon = [w.commonName[1] for w in warblers]
# Create index Dict for species names
speindex = indexmap(spenames)

# Observed coordinates range
coords_obs = (left = minimum(df.longitude), right = maximum(df.longitude),
              bottom = minimum(df.latitude), top = maximum(df.latitude))

## Get environmental data (with different training resolutions)
# WorldClim data
wc_vars = map(x -> worldclim(x, resolution = 10.0)[coords], [1,12]);
# Landcover data
lc_vars = map(x -> landcover(x, resolution = 10.0)[coords], 1:10);
# Training data with finer resolution
wc_vars_train = map(x -> worldclim(x, resolution = 5.0)[coords_obs], [1,12]);
lc_vars_train = map(x -> landcover(x, resolution = 5.0)[coords_obs], 1:10);

# Combine environmental data
env_vars = vcat(wc_vars, lc_vars)
env_vars_train = vcat(wc_vars_train, lc_vars_train)

## Get distribution for all species
# Runs only if create_distributions = true, as it can take a while. Loaded from file otherwise
# create_distributions = true # should distributions be computed (optional, loaded otherwise)
if (@isdefined create_distributions) && create_distributions == true
    @info "$(outcome) distributions to be created"
    # Select function to run given desired outcome
    if outcome == "raw"
        # Get raw distributions
        @time distributions = @showprogress pmap(x -> presence_absence(x, env_vars[1]), warblers)
        # @time distributions = @showprogress pmap(x -> presence_absence(x, env_vars_train[1], full_range = true, binary = false), warblers)
    elseif outcome == "bio"
        # Get sdm distributions (with different training resolutions)
        @time distributions = @showprogress pmap(x -> bioclim(x, env_vars, training_layers = env_vars_train), warblers);
    end
end

## Export distributions
# save_data = true # should data files be overwritten (optional)
if (@isdefined save_data) && save_data == true
    ## Export to JLD2
    # Export data to JLD2
    @info "Exporting data to JLD2 file ($(outcome) distributions data)"
    @save joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions spenames specommon speindex
    # Export to ZIP archive
    @info "Exporting data from JLD2 to ZIP archive ($(outcome) distributions data)"
    _zip_jld2(joinpath("data", "jld2", "$(outcome)-distributions.zip"),
              joinpath("data", "jld2", "$(outcome)-distributions.jld2"))
    # Make sure JLD2 timestamp is more recent than ZIP archive
    touch(joinpath("data", "jld2", "$(outcome)-distributions.jld2"))

    ## Export to CSV as Y matrix
    @info "Exporting data to CSV as Y matrix ($(outcome) distributions data)"
    # Get Y matrix
    Y = calculate_Y(distributions)
    inds_obs = _indsobs(Y)
    Yobs = _Yobs(Y, inds_obs)
    # Convert to dataframe
    spe_df = DataFrame(Yobs)
    rename!(spe_df, Symbol.("sp", 1:ncol(spe_df)))
    insertcols!(spe_df, 1, site = inds_obs)
    # Export to CSV
    CSV.write(joinpath("data", "proc", "distributions_spe_full.csv"), spe_df, delim="\t")
else
    # Load data
    @info "Data imported from file ($(outcome) distributions data)"
    @load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions spenames specommon speindex
end

## Count sites with presence per species
pres_counts = [length(filter(x -> x > 0.0, species.grid)) for species in distributions]
sort(pres_counts)

## Plot result
# Species 1
sp1 = "Setophaga_townsendi"
map_sp1 = plotSDM(distributions[speindex[sp1]], c=:BuPu,
                  title = "$(sp1) distribution ($(outcome))",
                  colorbar=:none, dpi=300)
scatter!(map_sp1, [NaN], label="Occurrence", color=:purple, markershape=:rect, markersize=2,
                        legend=:bottomright, legendfontsize=5)
# Species 2
sp2 = "Setophaga_petechia"
map_sp2 = plotSDM(distributions[speindex[sp2]], c=:BuPu,
                  title = "$(sp2) distribution ($(outcome))",
                  colorbar=:none, dpi=300)
scatter!(map_sp2, [NaN], label="Occurrence", color=:purple, markershape=:rect, markersize=2,
                        legend=:bottomright, legendfontsize=5)

## Export figures
# save_figures = true # should figures be overwritten (optional)
if (@isdefined save_figures) && save_figures == true
    @info "Figures saved ($(outcome) distributions)"
    savefig(map_sp1, joinpath("fig", outcome, "01_$(outcome)_sp-$(sp1).png")
    savefig(map_sp2, joinpath("fig", outcome, "01_$(outcome)_sp-$(sp2).png")
else
    @info "Figures not saved ($(outcome) distributions)"
end
