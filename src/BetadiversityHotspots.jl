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

include(joinpath("lib", "analysis.jl"))
export calculate_Y, _indsobs, _indsnotobs, _Yobs, _Ytransf,
       calculate_richness, calculate_lcbd

include(joinpath("lib", "beta-div.jl"))
export BD

include(joinpath("lib", "bioclim.jl"))
export bioclim

include(joinpath("lib", "csvdata.jl"))
export prepare_ebd_data, prepare_ebd_data!

include(joinpath("lib", "landcover.jl"))
export landcover

include(joinpath("lib", "overloads.jl"))
export clip

include(joinpath("lib", "plotSDM.jl"))
export plotSDM, plotSDM2, plot

include(joinpath("lib", "presence-absence.jl"))
export presence_absence

include(joinpath("lib", "quantiles.jl"))
export quantiles, quantiles!

include(joinpath("lib", "shapefiles.jl"))
# export download_shapefile, worldshape, clip, isin

include(joinpath("lib", "zipfile.jl"))
export verify_jld2_data, _unzip_jld2, _zip_jld2

end
