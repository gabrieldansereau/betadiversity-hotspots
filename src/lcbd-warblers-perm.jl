using Distributed
using JLD2
addprocs(9)

@everywhere using Random

@time @everywhere include("src/required.jl")

## Load functions to compute BD statistics
@everywhere include("src/lib/beta-div.jl")

# Load predictions for all species
@load "../data/predictions-can.jld2" predictions

## Create matrix Y (site-by-species community data table)
begin
    # Create Y -> site-by-species community data table
    Y = zeros(Int64, (prod(size(predictions[1])),length(predictions)))
    # Fill Y with community predictions
    @progress for gc in eachindex(predictions[1].grid) # loop for all sites
        # Group predictions for all species in site
        R = map(x -> x.grid[gc], predictions)
        # Fill Y with binary values -> 1 if species prediction for site != NaN, 0 if == NaN
        global Y[gc,:] = .!isnan.(R)
    end
end

#### 4 Options to calculate BD

## Option 1: Normal BD calculation (original one, same as lcbd-warblers.jl)
# Compute BD statistics
@time resBDnorm = BD(Y)

## Option 2: Calculate LCBD only for sites with predictions
# Get index of sites with predictions
sites_pred = map(x -> any(x .> 0), eachrow(Y))
inds_pred = findall(sites_pred)

# Select sites with predictions only
Ypred = Y[inds_pred,:]
# Compute BD statistics
resBDpred = BD(Ypred)

## Options 3-4: Permutation tests on options 1-2
nperm = 99
# Compute BD Statistics
@time resBDperm = BDperm(Y, nperm = nperm, distributed = true)
@time resBDpredperm = BDperm(Ypred, nperm = nperm, distributed = true)

## Extract LCBD values
resBD = [resBDnorm, resBDpred, resBDperm, resBDpredperm]
LCBDsets = [res.LCBDi for res in resBD]
# Scale LCBDi values to maximum value
LCBDsets = [LCBDi./maximum(LCBDi) for LCBDi in LCBDsets]

## Arrange LCBD values as grid
# Create empty grids
t_lcbd = [fill(NaN, size(predictions[1])) for LCBDi in LCBDsets]
# Fill in grid (options 1-2)
t_lcbd[1][inds_pred] = LCBDsets[1][inds_pred]
t_lcbd[2][inds_pred] = LCBDsets[2]
# Get indices of sites with significant LCBDs
inds_signif = [falses(prod(size(predictions[1]))) for LCBDi in LCBDsets]
inds_signif[3] = Array{Bool}(resBD[3].pLCBD .<= 0.05)
inds_signif[4][inds_pred] = Array{Bool}(resBD[4].pLCBD .<= 0.05)
# Get indices of sites with significant LCBDs & predictions
inds_signifpred = [intersect(findall(inds), inds_pred) for inds in inds_signif]
# Fill in grid (options 3-4)
[t_lcbd[i][inds_pred] .= 0.0 for i in 3:4]
[t_lcbd[i][inds_signifpred[i]] .= 1.0 for i in 3:4]

# Create SDMLayer with LCBD values
LCBD = SDMLayer.(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)

## Plot result
sdm_plot1 = plotSDM(LCBD[1], type="lcbd")
title!(sdm_plot1, "Option 1 - LCBD values with NaN sites")
sdm_plot2 = plotSDM(LCBD[2], type="lcbd")
title!(sdm_plot2, "Option 2 - LCBD values without NaN sites")
sdm_plot3 = plotSDM(LCBD[3], type="lcbd")
title!(sdm_plot3, "Option 3 - Significant LCBDs with NaN sites")
sdm_plot4 = plotSDM(LCBD[4], type="lcbd")
title!(sdm_plot4, "Option 4 - Significant LCBDs without NaN sites")
plots = plot(sdm_plot1, sdm_plot2, sdm_plot3, sdm_plot4, size=(1000,600))

## Save result
savefig(plots, "fig/warblers/lcbd-map-can-4options.pdf")
