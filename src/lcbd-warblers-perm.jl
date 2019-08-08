using Distributed
using JLD2

@time @everywhere include("src/required.jl")

## Get & prepare data
@time @everywhere begin
    # Load data from CSV files
    df = CSV.read("../data/warblers_can.csv", header=true, delim="\t")
    # Prepare data (select columns, arrange values)
    df = prepare_csvdata(df)
    # Separate species
    taxa_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-136.0, -58.0)
    lat_range = (40.5, 56.0)
end

# Load predictions for all species
@load "../data/predictions-can.jld2" predictions

## Create matrix Y (site-by-species community data table)
begin
    # Create Y -> site-by-species community data table
    Y = zeros(Int64, (prod(size(predictions[1])),length(taxa_occ)))
    # Fill Y with community predictions
    @progress for gc in eachindex(predictions[1].grid) # loop for all sites
        # Group predictions for all species in site
        R = map(x -> x.grid[gc], predictions)
        # Fill Y with binary values -> 1 if species prediction for site != NaN, 0 if == NaN
        global Y[gc,:] = .!isnan.(R)
    end

    # Get index of sites with predictions
    inds_pred = map(x -> any(x .> 0), eachrow(Y))
    # Select sites with predictions only
    Ypred = Y[inds_pred,:]
end

## Compute beta diversity statistics
# Load functions
@everywhere include("src/lib/beta-div.jl")
# Compute BD statistics
resBD = BD(Y)
# Extract LCBD values
LCBDi = resBD.LCBDi

## Arrange LCBD values as grid
# Create empty grid
t_lcbd = zeros(Float64, size(predictions[1]))
# Scale LCBDi values to maximum value
LCBDi = LCBDi./maximum(LCBDi)
t_lcbd[inds_pred] = LCBDi
t_lcbd[.!inds_pred] .= NaN
# Fill-in grid
#=
for i in eachindex(t_lcbd)
    # Add LCBD values only if prediction != NaN
    t_lcbd[i] = any(Y[i,:] .> 0) && p_LCBD[i] < 0.05 ? 1 : NaN # ?: is ternary operator, ~if-else
end
=#
# res_perm = BDperm(Ypred)
inds_signif = falses(prod(size(predictions[1])))
inds_signif[inds_pred] = Array{Bool}(res_perm.pLCBD .<= 0.05)
t_lcbd = zeros(Float64, size(predictions[1]))
t_lcbd[inds_signif] .= 1.0
t_lcbd[.!inds_signif] .= NaN
# Create SDMLayer with LCBD values
LCBD = SDMLayer(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot result
plotSDM(LCBD, type="lcbd")

## Save result
savefig(sdm_plot, "fig/warblers/lcbd-map-can-perm.pdf")
