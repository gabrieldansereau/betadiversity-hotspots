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

if nprocs() == 1
    include("BetadiversityHotspots.jl")
else
    @everywhere include(joinpath("src", "BetadiversityHotspots.jl"))
end
using .BetadiversityHotspots

verify_jld2_data(joinpath("data", "jld2"))