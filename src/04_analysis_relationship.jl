import Pkg; Pkg.activate(".")
using Distributed
@time @everywhere include(joinpath("src", "required.jl"))

## Conditional arguments
outcome = "raw"
# save_figures = true # should figures be overwritten (optional)

# Load distributions
@load joinpath("data", "jld2", "$(outcome)-distributions.jld2") distributions

# Calculate values
Y = calculate_Y(distributions)
richness = calculate_richness(Y, distributions[1])
lcbd = calculate_lcbd(Y, distributions[1])

## Plot relationship as histogram2d
rel2d = histogram2d(richness, lcbd, c = :viridis, bins = 40, title = "Relationship",
                    xlabel = "Richness", ylabel = "LCBD", colorbar_title = "Number of sites",
                    xlim = (1, 45), ylim = (0.0, 1.0), dpi = 150)

if (@isdefined save_figures) && save_figures == true
    savefig(rel2d, joinpath("fig", outcome, "06_$(outcome)_relationship2d-transf.png"))
end
