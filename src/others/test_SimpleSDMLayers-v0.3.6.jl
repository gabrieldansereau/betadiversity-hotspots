import Pkg
Pkg.activate(".")

# Pkg.develop("~/github/SimpleSDMLayers.jl/")

using SimpleSDMLayers
using JLD2
using Plots

@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
distrib = distributions[1]

temperature = worldclim(1)
temperature[distrib]
SimpleSDMLayers._layers_are_compatible(temperature, distrib)
test = distrib[distributions[2]]
isequal(distrib, test)
isequal(distrib.grid, test.grid) # true??
isequal(distrib.left, test.left)
isequal(distrib.right, test.right)
isequal(distrib.bottom, test.bottom)
isequal(distrib.top, test.top)
isequal(size(distrib), size(test))
# All true??
# whatever...

coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
temperature = temperature[coords_qc]
distrib = distrib[coords_qc]

## Testing new functions

# similar
similar(temperature)
similar(Float64, temperature)

# hcat/vcat
l1 = worldclim(1, left=0.0, right=10.0, bottom=0.0, top=10.0)
l2 = worldclim(1, left=0.0, right=10.0, bottom=10.0, top=20.0)
l3 = worldclim(1, left=10.0, right=20.0, bottom=0.0, top=10.0)
l4 = worldclim(1, left=10.0, right=20.0, bottom=10.0, top=20.0)

ml1 = hcat(l1, l3)
vl1 = vcat(l1, l2)
ml2 = hcat(l2, l4)
vl2 = vcat(l3, l4)

# collect
collect(distrib)
collect(temperature)
isequal(collect(temperature), filter(!isnothing, temperature.grid))

# isnothing
isnothing(distrib)

# mask
SimpleSDMLayers._inner_type(temperature)
mask(distrib, temperature)
plot(temperature, c = :viridis)
plot(distrib, c = :viridis)
plot(mask(distrib, temperature), c = :viridis)

# plot recipes
heatmap(temperature)
scatter(distrib, temperature)
