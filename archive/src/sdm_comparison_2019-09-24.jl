using Distributed
using JLD2
@time include("../required.jl")

# Precompile plot function
@time plot()

# Load results to compare
@time include("04_sdm_community.jl")
@time include("05_sdm_lcbd.jl")
@time include("06_sdm_richness.jl")

# Combine plots
comparison_plot = plot(diversity_plot, lcbd_plot,
                       richness_plot, richness_zscores_plot,
                       size=(1000,700))

# Save results
#=
savefig(comparison_plot, "fig/sdm/sdm-comparison.pdf")
=#
