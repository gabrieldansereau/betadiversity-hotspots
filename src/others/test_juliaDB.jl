using JuliaDB
using IndexedTables
using Dates

@time df = loadtable("data/raw/ebd_warblers_cut.csv", delim='\t',
                        nastrings=[""],
                        type_detect_rows=1000, spacedelim=false)

# Fix names case & spacing
oldnames = collect(colnames(df))
newnames = oldnames .|>
    string .|>
    titlecase .|>
    lowercasefirst .|>
    x -> replace(x, "_" => "") .|>
    Symbol
replace!(newnames, :scientificName => :species)
df = rename(df, oldnames .=> newnames)
# Separate year-month-day
years = year.(select(df, :observationDate))
df = transform(df, :year => years)
# Remove entries with missing data
IndexedTables.dropmissing(df, :species)
dropmissing(df)
