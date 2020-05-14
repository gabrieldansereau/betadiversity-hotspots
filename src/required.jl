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

@everywhere include(joinpath("src", "BetadiversityHotspots.jl"))
using .BetadiversityHotspots
