## Investigation in progress !!!

# distributions_spe_full.csv has changed when allowing type Nothing
# but it shouldn't have!!
# why is that?

using CSV
using DelimitedFiles

path_new = joinpath("data", "proc", "distributions_spe_full.csv")
path_old = joinpath("data", "proc", "distributions_spe_full_backup.csv")

# Checking the CSV files
Ynew = readdlm(path_new, header = true)[1]
Yold = readdlm(path_old, header = true)[1]
spe_new = CSV.File(path_new, header = true, delim = "\t") |> DataFrame!
spe_old = CSV.File(path_old, header = true, delim = "\t") |> DataFrame!
spe_new
spe_old
# different number of lines and columns ...

## Investigating different number of columns

names(spe_new)
names(spe_old)
symdiff(names(spe_new), names(spe_old))
# ok so sp63 is the difference... where does it come from??

# Checking the JLD2 files
outcome = "raw"
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions
distributions_new = distributions
@load joinpath("data", "jld2", "$(outcome)-distributions_nan.jld2") distributions
distributions_old = distributions
# 63 species!!
[any(x -> !isnan(x), d.grid) for d in distributions_old]
# sp60 only has NaN!

# The filter() call was moved earlier in the workflow in commit 174e01 "something broken with filter"
# It is now performed before calling `presence_absence()`, so sp60 is removed beforehand

# Test internals of presence_absence()
testdf = DataFrame(species = ["test"], longitude = [255.0], latitude = [255.0])
coords = (left = -145.0, right = -50.0, bottom = 20.0, top = 75.0)
filter!(x -> coords.left < x.longitude < coords.right, testdf)
unique_sites = unique(testdf, [:longitude, :latitude])
# Works with empty DataFrame, so the presence_absence() call didn't fail,
# and returned an empty layer for sp60
# (which was removed before the predictions in R I guess?)

# Solved!

## Investigating different number of sites 

# Get different sites
diffsites = symdiff(spe_new.site, spe_old.site)
# 67 sites different

# Plot'em
temp = worldclim(1)[coords]
difflayer = similar(temp)
replace!(x -> isnothing(x) ? x : nothing, difflayer.grid)
difflayer.grid[diffsites] .= 1.0
plotSDM2(difflayer)

# Get coordinates
spadf = CSV.File(joinpath("data", "proc", "distributions_spa_full.csv"), header = true) |> DataFrame!
show(spadf[diffsites, :], allrows = true)
# 60 sites are just outside the south border (lat = 20.0)
# 6 sites aure just outside the west border (lon = -145.0)
# 1 site is just weird, guess it's a precision error

# Observations for these border sites should be included, I'll modify
# distributions.jl accordingly

# Solved!

