using Distributed
using JLD2
@time include("required.jl")

# Precompile plot function
@time plot()

# Load results to compare
@time include("community-warblers.jl")
@time include("lcbd-warblers.jl")
@time include("richness-warblers.jl")

# Combine plots
comparison_plot = plot(diversity_plot, lcbd_plot,
                       richness_plot, richness_zscores_plot,
                       size=(1000,700))

# Save results
savefig(comparison_plot, "fig/comparison-am-larger2.pdf")
