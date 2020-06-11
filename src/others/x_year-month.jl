import Pkg
Pkg.activate(".")
using Distributed
@time include("required.jl")

# Load data from CSV files (from file cut with terminal)
@time df = CSV.read(joinpath("data", "raw", "ebd_warblers_cut.csv"), header=true, delim="\t")

# Prepare data (arrange values & columns)
df = prepare_ebd_data(df)

# Remove duplicates (BUT not ok for counts, so-so for dates, see group-observation.jl in src/test/)
df_nogroups = filter(x -> ismissing(x[:groupIdentifier]), df)
df_groups = dropmissing(df, :groupIdentifier)
unique!(df_groups, [:species, :groupIdentifier])
append!(df_nogroups, df_groups)
newdf = df_nogroups

ycount = by(newdf, :year,  n = :year  => length, sort = true)
show(ycount, allrows = true)
bar(ycount.year, ycount.n)
filter(x -> x.year > 1980, ycount) |> x -> bar(x.year, x.n)
filter(x -> x.year > 1980, ycount) |> x -> bar(x.year, x.n, yaxis = :log)
filter(x -> x.year < 2000, ycount) |> x -> bar(x.year, x.n)
filter(x -> 1980 <= x.year <= 2000, ycount) |> x -> bar(x.year, x.n)

by(newdf, [:year, :protocolType],  n = :year  => length, sort = true)

mcount = by(newdf, :month, n = :month => length, sort = true)
bar(mcount.month, mcount.n)
