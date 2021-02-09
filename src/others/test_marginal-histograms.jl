using SimpleSDMLayers
using DataFrames
using RDatasets
using StatsPlots

iris = dataset("datasets","iris")
@df iris marginalhist(:PetalLength, :PetalWidth, bins = 40)

coords_qc = (left = -80.0, right = -55.0, bottom = 45.0, top = 63.0)
temperature = worldclim(1)[coords_qc]
temperature_df = DataFrame(temperature)

@df temperature_df marginalhist(:longitude, :latitude, bins = 40)

plot(temperature, c = :lightgrey)
@df temperature_df marginalhist!(:longitude, :latitude, bins = 40)

p = @df temperature_df marginalhist(:longitude, :latitude, bins = 40)
plot!(temperature, subplot = 2, aspectratio = :none)
plot!(aspectratio = [0.002 92/60 :auto])

savefig("testplot.png")

using JLD2
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions
distrib = distributions[1][coords_qc]

density(distrib)
density(temperature)


