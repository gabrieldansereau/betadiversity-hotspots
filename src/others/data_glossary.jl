import Pkg
Pkg.activate(".")
using Distributed
@time include(joinpath("..", "required.jl"))

## Create glossary dataframe

# Landcover data
lcnames = ["bare", "crops", "grass", "moss", "shrub", 
           "snow", "tree", "urban", "water-permanent", "water-seasonal"]
lcdf = DataFrame(variable = string.("lc", eachindex(lcnames)),
                 type = "landcover",
                 value = lcnames)

# Climate data
wcnames = ["Annual Mean Temperature"
           "Mean Diurnal Range (Mean of monthly (max temp - min temp))"
           "Isothermality (BIO2/BIO7) (* 100)"
           "Temperature Seasonality (standard deviation *100)"
           "Max Temperature of Warmest Month"
           "Min Temperature of Coldest Month"
           "Temperature Annual Range (BIO5-BIO6)"
           "Mean Temperature of Wettest Quarter"
           "Mean Temperature of Driest Quarter"
           "Mean Temperature of Warmest Quarter"
           "Mean Temperature of Coldest Quarter"
           "Annual Precipitation"
           "Precipitation of Wettest Month"
           "Precipitation of Driest Month"
           "Precipitation Seasonality (Coefficient of Variation)"
           "Precipitation of Wettest Quarter"
           "Precipitation of Driest Quarter"
           "Precipitation of Warmest Quarter"
           "Precipitation of Coldest Quarter"]
wcdf = DataFrame(variable = string.("wc", eachindex(wcnames)),
                 type = "climate",
                 value = wcnames)

# Species data
@load joinpath("data", "jld2", "raw-distributions.jld2") spenames
spenames
spdf = DataFrame(variable = string.("sp", eachindex(spenames)),
                 type = "species",
                 value = spenames)

# Create glossary
glossdf = vcat(lcdf, wcdf, spdf)

# Export to CSV
CSV.write(joinpath("data", "proc", "glossary.csv"), glossdf)