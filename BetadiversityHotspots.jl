module BetadiversityHotspots

using JLD2
using Plots
using Shapefile
using GBIF
using StatsBase
using Statistics
using DataFrames
using CSV
using DelimitedFiles
using Dates
using SimpleSDMLayers
using Random
using ProgressMeter
using ArchGDAL
using RCall
using ZipFile

include(joinpath("src", "lib", "analysis.jl"))
export calculate_Y, _indsobs, _indsnotobs, _Yobs, _Ytransf,
       calculate_richness, calculate_lcbd

include(joinpath("src", "lib", "beta-div.jl"))
export BD

include(joinpath("src", "lib", "bioclim.jl"))
export _bcscore, _bioclim_layer, bioclim

include(joinpath("src", "lib", "csvdata.jl"))
export prepare_gbif_data, prepare_ebd_data

include(joinpath("src", "lib", "landcover.jl"))
export landcover

include(joinpath("src", "lib", "overloads.jl"))
export getindex, longitudes, latitudes, clip, minimum

include(joinpath("src", "lib", "plotSDM.jl"))
export plotSDM, plot

include(joinpath("src", "lib", "presence-absence.jl"))
export presence_absence

include(joinpath("src", "lib", "quantiles.jl"))
export quantiles, quantiles!

include(joinpath("src", "lib", "shapefiles.jl"))
# export download_shapefile, worldshape, clip, isin

include(joinpath("src", "lib", "zipfile.jl"))
export verify_jld2_data, _unzip_jld2, _zip_jld2


end
