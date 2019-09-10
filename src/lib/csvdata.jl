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
