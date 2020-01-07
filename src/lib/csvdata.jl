function prepare_gbif_data(df::DataFrame)
    # Rename coordinate columns names
    rename!(df, :decimalLatitude => :latitude)
    rename!(df, :decimalLongitude => :longitude)
    # Select subset with specific columns
    select!(df, [:species, :year, :latitude, :longitude])
    # Remove entries with missing data
    df = dropmissing(df, :year)
    df = dropmissing(df, :species)
    # Replace spaces by underscores in species names
    df.species .= replace.(df.species, " " .=> "_")
    return df
end

function prepare_ebd_data(df::DataFrame)
    # Fix names case & spacing
    newnames = names(df) .|>
        string .|>
        titlecase .|>
        lowercasefirst .|>
        x -> replace(x, " " => "") .|>
        Symbol
    names!(df, newnames)
    # Rename species column
    rename!(df, :scientificName => :species)
    # Separate year-month-day
    df.year = year.(df.observationDate)
    # Remove entries with missing data
    df = dropmissing(df, :species)
    # Replace spaces by underscores in species names
    df.species .= replace.(df.species, " " .=> "_")
    # Remove not approved observations (exotic species, unvetted data)
    filter!(obs -> obs[:approved] .== 1, df)
    return df
end

function load_landcover(lon_range, lat_range;path::AbstractString="assets/landcover/")
    # Load variable names & paths
    landcover_files = readdir(path)
    landcover_vars = [split(lc, r"_")[2] for lc in landcover_files]
    landcover_paths = string.(path, landcover_files)

    # Load data
    landcover_raw = readdlm.(landcover_paths)
    # Reverse rows
    landcover_mat = [lc[end:-1:1, :] for lc in landcover_raw]

    # Replace 255 (default no data values) by NaN
    [replace!(lc, 255 => NaN) for lc in landcover_mat]
    # Add additionnal line of NaNs down South (coordinates prime did not exactly match)
    landcover_mat = [vcat(fill(NaN, (1, size(lc, 2))), lc) for lc in landcover_mat]

    # Convert to SDMLayers, cut to coordinates range
    landcover_layers = [SimpleSDMPredictor(lc, -160.0, -40.0, 20.0, 80.0)[lon_range, lat_range] for lc in landcover_mat]
end
