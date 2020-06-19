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
wcnames = ["Annual Mean Temperature (°C)"
           "Mean Diurnal Range (Mean of monthly (max temp - min temp)) (°C)"
           "Isothermality (BIO2/BIO7) (* 100) (°C)"
           "Temperature Seasonality (standard deviation *100) (°C)"
           "Max Temperature of Warmest Month (°C)"
           "Min Temperature of Coldest Month (°C)"
           "Temperature Annual Range (BIO5-BIO6) (°C)"
           "Mean Temperature of Wettest Quarter (°C)"
           "Mean Temperature of Driest Quarter (°C)"
           "Mean Temperature of Warmest Quarter (°C)"
           "Mean Temperature of Coldest Quarter (°C)"
           "Annual Precipitation (mm)"
           "Precipitation of Wettest Month (mm)"
           "Precipitation of Driest Month (mm)"
           "Precipitation Seasonality (Coefficient of Variation) (mm)"
           "Precipitation of Wettest Quarter (mm)"
           "Precipitation of Driest Quarter (mm)"
           "Precipitation of Warmest Quarter (mm)"
           "Precipitation of Coldest Quarter (mm)"]
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