using Distributed
using JLD2
using Plots
using Plots.PlotMeasures
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
using Formatting
using StatsPlots

include(joinpath("lib", "analysis.jl"))
include(joinpath("lib", "beta-div.jl"))
include(joinpath("lib", "bioclim.jl"))
include(joinpath("lib", "csvdata.jl"))
include(joinpath("lib", "landcover.jl"))
include(joinpath("lib", "plotSDM.jl"))
include(joinpath("lib", "presence-absence.jl"))
include(joinpath("lib", "quantiles.jl"))
include(joinpath("lib", "shapefiles.jl"))
include(joinpath("lib", "zipfile.jl"))

# if nprocs() == 1
#     include("BetadiversityHotspots.jl")
# else
#     @everywhere include(joinpath("src", "BetadiversityHotspots.jl"))
# end
# using .BetadiversityHotspots

verify_jld2_data(joinpath("data", "jld2"))