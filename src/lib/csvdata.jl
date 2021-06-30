## Data preparation functions

# eBird data preparation (from DataFrame)

"""
See `prepare_ebd_data!` for details.
"""
function prepare_ebd_data(df::DataFrame)
    # Fix names case & spacing
    newnames =
        names(df) .|>
        string .|>
        titlecase .|>
        lowercasefirst .|>
        x -> replace(x, " " => "") .|> Symbol
    df = rename(df, newnames)
    # Rename species column
    df = rename(df, :scientificName => :species)
    # Separate year-month-day
    df.year = year.(df.observationDate)
    df.month = month.(df.observationDate)
    # Remove entries with missing data
    df = dropmissing(df, :species)
    # Replace spaces by underscores in species names
    df.species .= replace.(df.species, " " .=> "_")
    # Remove not approved observations (exotic species, unvetted data)
    df = filter(:approved => ==(1), df)
    return df
end

"""
    prepare_ebd_data!(df::DataFrame)
    prepare_ebd_data(df::DataFrame)

Prepares the eBird data for the analyses from this project. Currently, this
fixes the column and species names, adds a year and a month column, removes
entries with missing data, and removed unapproved observations (exotic species,
unvetted data).

`prepare_ebd_data(df::DataFrame)` copies the DataFrame, while
`prepare_ebd_data!(df::DataFrame)` operates in-place. The in-place version is
preferred when working on the very large eBird DataFrame.
"""
function prepare_ebd_data!(df::DataFrame)
    # Fix names case & spacing
    newnames =
        names(df) .|>
        string .|>
        titlecase .|>
        lowercasefirst .|>
        x -> replace(x, " " => "") .|> Symbol
    rename!(df, newnames)
    # Rename species column
    rename!(df, :scientificName => :species)
    # Separate year-month-day
    df.year = year.(df.observationDate)
    df.month = month.(df.observationDate)
    # Remove entries with missing data
    dropmissing!(df, :species)
    # Replace spaces by underscores in species names
    df.species .= replace.(df.species, " " .=> "_")
    # Remove not approved observations (exotic species, unvetted data)
    filter!(:approved => ==(1), df)
    return df
end
