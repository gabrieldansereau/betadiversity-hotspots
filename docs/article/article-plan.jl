import Pkg
Pkg.activate(".")
using Distributed
@time include(joinpath("..", "..", "src", "required.jl"))

outcome = "bart"

# Run analyses
@time include(joinpath("..", "..", "src", "04_analysis.jl"))

# Remove figures titles
plot!(richness_plot, title = "")
plot!(lcbdtr_plot, title = "", clim = (-Inf, Inf), cbartitle = "Relative LCBD value")
plot!(rel2d_plot, title = "")

# Export figures
savefig(plot(richness_plot, dpi = 150), joinpath("docs", "article", "fig", "richness.png"))
savefig(plot(lcbdtr_plot, dpi = 150), joinpath("docs", "article", "fig", "lcbd.png"))
savefig(plot(rel2d_plot, dpi = 150), joinpath("docs", "article", "fig", "relationship.png"))

