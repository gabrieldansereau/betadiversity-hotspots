include("../required.jl")

outcome = "bart"
outcome = "raw"

files = [
    "fig/$(outcome)/04-3_$(outcome)_lcbd-transf.png"                 "fig/$(outcome)/04_$(outcome)_lcbd.png";
    "fig/$(outcome)/04-4_$(outcome)_relationship2d-transf.png"       "fig/$(outcome)/04_$(outcome)_relationship.png";
    "fig/$(outcome)/04-2_$(outcome)_richness.png"                    "fig/$(outcome)/04_$(outcome)_richness.png";
    "fig/$(outcome)/05-2_$(outcome)_subareas_3scales.png"            "fig/$(outcome)/05_$(outcome)_extents.png";
    "fig/$(outcome)/05-4_$(outcome)_subareas_medians.png"            "fig/$(outcome)/05_$(outcome)_medians.png";
    "fig/$(outcome)/05-3_$(outcome)_subareas.gif"                    "fig/$(outcome)/05_$(outcome)_scaling.gif";
    "fig/$(outcome)/05-1_$(outcome)_subareas_combined.png"           "fig/$(outcome)/05_$(outcome)_subareas.png";
    "fig/$(outcome)/07_$(outcome)_combined-maps.png"                 "fig/$(outcome)/09_$(outcome)_combined.png";
    "fig/$(outcome)/07_$(outcome)_comparison-combined.png"           "fig/$(outcome)/09_$(outcome)_difference.png";
    "fig/$(outcome)/07_$(outcome)_residuals-combined.png"            "fig/$(outcome)/09_$(outcome)_residuals.png";
    "fig/$(outcome)/08_$(outcome)_rare-species_eusrr_thresholds.png" "fig/$(outcome)/06_$(outcome)_eusrr.png";
    "fig/$(outcome)/08_$(outcome)_rare-species_ascending_plots.png"  "fig/$(outcome)/06_$(outcome)_rare-species.png";
    "fig/$(outcome)/08_$(outcome)_rare-species_thresholds.png"       "fig/$(outcome)/06_$(outcome)_thresholds.png";
    "fig/$(outcome)/06-0_$(outcome)_moving-windows_full.png"          "archive/fig/$(outcome)/06-0_$(outcome)_moving-windows_full.png";
    "fig/$(outcome)/06-1_$(outcome)_moving-windows_subareas.png"      "archive/fig/$(outcome)/06-1_$(outcome)_moving-windows_subareas.png";
    "fig/$(outcome)/06-2_$(outcome)_moving-windows_3scales.png"       "archive/fig/$(outcome)/06-2_$(outcome)_moving-windows_3scales.png";
    "fig/$(outcome)/06-3_$(outcome)_moving-windows.gif"               "archive/fig/$(outcome)/06-3_$(outcome)_moving-windows.gif";
    "fig/$(outcome)/08_$(outcome)_rare-species_spatial_subareas.png"  "archive/fig/$(outcome)/08_$(outcome)_rare-species_spatial_subareas.png";
]

oldfiles = files[:, 1]
newfiles = files[:, 2]

isfile.(oldfiles)
filter(!isfile, oldfiles)

for (old, new) in zip(oldfiles, newfiles)
    if isfile(old)
        mv(old, new)
    end
end

isfile.(newfiles)
filter(!isfile, newfiles)

oldfiles = replace.(oldfiles, "fig/$(outcome)/" => "")
newfiles = replace.(newfiles, "archive/fig/$(outcome)/" => "")
newfiles = replace.(newfiles, "fig/$(outcome)/" => "")

file = "src/04_full-extent.jl"
(tmppath, tmpio) = mktemp()
cp(file, tmppath; force=true)
for (old, new) in zip(oldfiles[1:3], newfiles[1:3])
    open(tmppath) do io
        for line in eachline(io, keep=true) # keep so the new line isn't chomped
            line = replace(line, old => new)
            write(tmpio, line)
        end
    end
    close(tmpio)
end
mv(tmppath, file, force=true)
