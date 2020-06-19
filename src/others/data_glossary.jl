import Pkg
Pkg.activate(".")
using Distributed
@time include(joinpath("..", "required.jl"))

## Create glossary dataframe

# Landcover data
lcgloss = ["lc1"    "bare"               "Bare land cover (%)"
           "lc2"    "crops"              "Crops land cover (%)"
           "lc3"    "grass"              "Grass land cover (%)"
           "lc4"    "moss"               "Moss land cover (%)"
           "lc5"    "shrub"              "Shrub land cover (%)"
           "lc6"    "snow"               "Snow land cover (%)"
           "lc7"    "tree"               "Tree land cover (%)"
           "lc8"    "urban"              "Urban land cover (%)"
           "lc9"    "water_permanent"    "Permanent water land cover (%)"
           "lc10"   "water_seasonal"     "Seasonal water land cover (%)"]
lcdf = DataFrame(variable = lcgloss[:,1],
                 type = "landcover",
                 full_name = lcgloss[:,2],
                 description = lcgloss[:,3])
# Climate data
wcgloss = ["wc1"     "temperature_annual_mean"             "Annual Mean Temperature (°C)"
           "wc2"     "temperature_diurnal_range"           "Mean Diurnal Range (Mean of monthly (max temp - min temp)) (°C)"
           "wc3"     "temperature_isothermality"           "Isothermality (BIO2/BIO7) (* 100) (°C)"
           "wc4"     "temperature_seasonality"             "Temperature Seasonality (standard deviation *100) (°C)"
           "wc5"     "temperature_max_warmest_month"       "Max Temperature of Warmest Month (°C)"
           "wc6"     "temperature_min_warmest_month"       "Min Temperature of Coldest Month (°C)"
           "wc7"     "temperature_annual_range"            "Temperature Annual Range (BIO5-BIO6) (°C)"
           "wc8"     "temperature_mean_wettest_quarter"    "Mean Temperature of Wettest Quarter (°C)"
           "wc9"     "temperature_mean_driest_quarter"     "Mean Temperature of Driest Quarter (°C)"
           "wc10"    "temperature_mean_warmest_quarter"    "Mean Temperature of Warmest Quarter (°C)"
           "wc11"    "temperature_mean_coldest_quarter"    "Mean Temperature of Coldest Quarter (°C)"
           "wc12"    "precipitation_annual_mean"           "Annual Precipitation (mm)"
           "wc13"    "precipitation_wettest_month"         "Precipitation of Wettest Month (mm)"
           "wc14"    "precipitation_driest_month"          "Precipitation of Driest Month (mm)"
           "wc15"    "precipitation_seasonality"           "Precipitation Seasonality (Coefficient of Variation) (mm)"
           "wc16"    "precipitation_wettest_quarter"       "Precipitation of Wettest Quarter (mm)"
           "wc17"    "precipitation_driest_quarter"        "Precipitation of Driest Quarter (mm)"
           "wc18"    "precipitation_warmest_quarter"       "Precipitation of Warmest Quarter (mm)"
           "wc19"    "precipitation_coldest_quarter"       "Precipitation of Warmest Quarter (mm)"]
wcdf = DataFrame(variable = wcgloss[:,1],
                 type = "climate",
                 full_name = wcgloss[:,2],
                 description = wcgloss[:,3])

# Species data
@load joinpath("data", "jld2", "raw-distributions.jld2") spenames
spenames
spdf = DataFrame(variable = string.("sp", eachindex(spenames)),
                 type = "species",
                 full_name = spenames,
                 description = replace.(spenames, "_" => " "))

# Create glossary
glossdf = vcat(lcdf, wcdf, spdf)

# Export to CSV
CSV.write(joinpath("data", "proc", "glossary.csv"), glossdf)