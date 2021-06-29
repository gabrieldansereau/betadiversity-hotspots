# Activate project
using Pkg: Pkg
Pkg.activate(".")
Pkg.instantiate()

# Load required packages (sorted alphabetically)
using CSV
using DataFrames
using Dates
using DelimitedFiles
using Distributed
using Formatting
using JLD2
using Plots
using Plots.PlotMeasures
using ProgressMeter
using Random
using RCall
using Shapefile
using SimpleSDMLayers
using Statistics
using StatsBase
using StatsPlots
using ZipFile

# Load custom functions
include(joinpath("lib", "analysis.jl"))
include(joinpath("lib", "beta-div.jl"))
include(joinpath("lib", "csvdata.jl"))
include(joinpath("lib", "landcover.jl"))
include(joinpath("lib", "plotSDM.jl"))
include(joinpath("lib", "presence-absence.jl"))
include(joinpath("lib", "shapefiles.jl"))

# Manual version-control for oversized files
include(joinpath("lib", "version-control.jl"))
verify_raw_files()
verify_proc_files()
# verify_raw_files(; touch_placeholders=true)
# verify_proc_files(; touch_placeholders=true)
