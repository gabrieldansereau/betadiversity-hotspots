using Pkg: Pkg
Pkg.activate(".")
using Distributed
@time @everywhere include("src/required.jl")

#### Richness on raw data (for the mask)

## Load presence-absence data for all species
@load "data/jld2/raw-pres-abs.jld2" pres_abs

## Load matrix Y
@load "data/jld2/raw-Y-matrices.jld2" Y Ypred Ytransf inds_pred inds_notpred

#### Species richness
## Get number of species per site
sums_raw = map(x -> Float64(sum(x)), eachrow(Y))
# Add NaN for non predicted sites
sums_raw[inds_notpred] .= NaN
# Reshape to grid format
sums_raw = reshape(sums_raw, size(pres_abs[1]))

#### Richness on sdm predictions

## Load predictions for all species
@load "data/jld2/sdm-predictions-sametrain.jld2" predictions

## Load matrix Y
@load "data/jld2/sdm-Y-matrices-sametrain.jld2" Y Ypred Ytransf inds_pred inds_notpred

#### Species richness
## Get number of species per site
sums = map(x -> Float64(sum(x)), eachrow(Y))
# Add NaN for non predicted sites
sums[inds_notpred] .= NaN
# Reshape to grid format
sums = reshape(sums, size(predictions[1]))

#### Apply mask
# Empty array to collect plots
mask_plots = []
# Select ranges to test
ranges = [1:10..., 15, 20, 25, 30]
# Apply mask using different ranges
@time @progress for r in ranges
    # Determine range
    ncells = r
    # Verify if any non-Nan value within the range
    verified = falses(size(sums))
    for i in 1:size(sums, 1), j in 1:size(sums, 2)
        # Define number of cells to check on each side
        ni_inf, ni_sup, nj_inf, nj_sup = repeat([ncells - 1], 4)
        # If-else to make things work with outer cells
        if i < ncells
            ni_inf = i - 1
        elseif i > size(sums, 1) - ncells
            ni_sup = size(sums, 1) - i - 1
        end
        if j < ncells
            nj_inf = j - 1
        elseif j > size(sums, 2) - ncells
            nj_sup = size(sums, 2) - j - 1
        end
        # Verify if any non-Nan value within the range
        verified[i, j] = any(
            !isnan, sums_raw[(i - ni_inf):(i + ni_sup), (j - nj_inf):(j + nj_sup)]
        )
    end
    # Create mask
    masked = .!verified
    sums_mask = copy(sums)
    # Change masked values to NaN
    sums_mask[masked] .= NaN
    # Create SimpleSDMLayer
    richness_mask = SimpleSDMResponse(
        sums_mask,
        predictions[1].left,
        predictions[1].right,
        predictions[1].bottom,
        predictions[1].top,
    )
    # Plot result
    richness_plot_mask = plotSDM(richness_mask; c=:viridis)
    title!(richness_plot_mask, "ncells = $(ncells)")
    # Export plot
    push!(mask_plots, richness_plot_mask)
end

# Create gif
anim = @animate for i in 1:length(mask_plots)
    plot(mask_plots[i])
end
gif(anim, "fig/sdm/sdm-richness-mask.gif"; fps=1)

## Save result
#=
savefig(richness_plot, "fig/sdm/03_sdm_richness.pdf")
savefig(richness_plot_mask, "fig/sdm/sdm-richness-mask.pdf")
=#
