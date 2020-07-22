if !(@isdefined BetadiversityHotspots)
    import Pkg; Pkg.activate(".")
    @time include(joinpath("..", "..", "src", "required.jl"))
end

# Run analyses
# outcome = "raw"
outcome = "bart"
@time include(joinpath("..", "..", "src", "04_analysis.jl"))

# Remove figures titles
plot!(richness_plot, title = "")
plot!(lcbdtr_plot, title = "", cbartitle = "Relative LCBD value")
plot!(rel2d_plot, title = "", 
      size = (900, 400),
      ratio = 40,
      tickfontsize = 8,
      # left_margin = 0.0Plots.PlotMeasures.mm, right_margin = 0.0Plots.PlotMeasures.mm,
      )

# Export figures
savefig(plot(richness_plot, dpi = 150), joinpath("docs", "article", "fig", "richness-$(outcome).png"))
savefig(plot(lcbdtr_plot, dpi = 150), joinpath("docs", "article", "fig", "lcbd-$(outcome).png"))
savefig(plot(rel2d_plot, dpi = 150), joinpath("docs", "article", "fig", "relationship-$(outcome).png"))
