function prepare_csvdata(csvdata::DataFrame)
    # Subset with specific columns
    df = csvdata[:, [:species, :year, :decimalLatitude, :decimalLongitude]]
    # Rename coordinate columns names
    rename!(df, :decimalLatitude => :latitude)
    rename!(df, :decimalLongitude => :longitude)
    # Remove entries with missing data
    dropmissing!(df, :year)
    dropmissing!(df, :species)
    # Replace spaces by underscores in species names
    df.species .= replace.(df.species, " " .=> "_")
    return df
end
