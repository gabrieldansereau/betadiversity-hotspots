## Speed test for CSV & DataFrames 
import Pkg;
Pkg.activate(".")

using CSV, DataFrames

## Test loading with depwarn
# Add "julia.additionalArgs": ["--depwarn=yes"] to settings.json
DataFrame!(CSV.File(joinpath("data", "proc", "bart_summaries.csv")));

## Very small DataFrame (6 kb)
@time DataFrame(CSV.File(joinpath("data", "proc", "bart_summaries.csv")));
# 8.417254 seconds (16.23 M allocations: 739.774 MiB, 5.39% gc time) # 1st ever
# 0.000371 seconds (311 allocations: 29.609 KiB) # 2nd
@time CSV.read(joinpath("data", "proc", "bart_summaries.csv"), DataFrame);
# 8.820570 seconds (16.24 M allocations: 740.719 MiB, 4.78% gc time) # 1st ever
# 0.000461 seconds (307 allocations: 26.359 KiB) # 2nd
# ---> equivalent

## Big DataFrame (1.9 GB)
@time DataFrame(CSV.File(joinpath("data", "proc", "ebd_warblers_prep.csv"), header=true, delim="\t"));
# 21.928007 seconds (90.90 M allocations: 7.392 GiB, 20.60% gc time) # 1st ever
# 19.391606 seconds (67.63 M allocations: 6.335 GiB, 84.98% gc time) # 2nd (doesn't get better)
@time CSV.read(joinpath("data", "proc", "ebd_warblers_prep.csv"), DataFrame, header=true, delim="\t");
# 19.678142 seconds (90.80 M allocations: 6.144 GiB, 14.62% gc time) # 1st ever
# 6.644478 seconds (67.63 M allocations: 5.091 GiB, 71.14% gc time) # 2nd (huge difference)
# 2nd way better on 2nd call

## Bunch of DataFrames, some relatively big
@time begin
    results  = DataFrame(CSV.File(joinpath("data", "proc", "bart_summaries.csv")));
    varimps  = DataFrame(CSV.File(joinpath("data", "proc", "bart_varimps.csv")));
    pred_df  = DataFrame(CSV.File(joinpath("data", "proc", "bart_predictions_prob.csv"), missingstrings = ["NA"]));
    lower_df = DataFrame(CSV.File(joinpath("data", "proc", "bart_predictions_lower.csv"), missingstrings = ["NA"]));
    upper_df = DataFrame(CSV.File(joinpath("data", "proc", "bart_predictions_upper.csv"), missingstrings = ["NA"]));
    pres_df  = DataFrame(CSV.File(joinpath("data", "proc", "bart_predictions_pres.csv"), missingstrings = ["NA"]));
end;
# 22.359386 seconds (33.94 M allocations: 2.324 GiB, 5.99% gc time) # 1st ever
# 22.450855 seconds (33.60 M allocations: 2.307 GiB, 5.05% gc time) # 1st ever without assigning
# 0.941610 seconds (105.62 k allocations: 852.245 MiB, 14.71% gc time) # 2nd
@time begin
    results  = CSV.read(joinpath("data", "proc", "bart_summaries.csv"), DataFrame)
    varimps  = CSV.read(joinpath("data", "proc", "bart_varimps.csv"), DataFrame)
    pred_df  = CSV.read(joinpath("data", "proc", "bart_predictions_prob.csv"), DataFrame, missingstrings = ["NA"])
    lower_df = CSV.read(joinpath("data", "proc", "bart_predictions_lower.csv"), DataFrame, missingstrings = ["NA"])
    upper_df = CSV.read(joinpath("data", "proc", "bart_predictions_upper.csv"), DataFrame, missingstrings = ["NA"])
    pres_df  = CSV.read(joinpath("data", "proc", "bart_predictions_pres.csv"), DataFrame, missingstrings = ["NA"])
end;
# 22.228993 seconds (33.98 M allocations: 1.977 GiB, 4.81% gc time) # 1st ever
# 21.759907 seconds (34.13 M allocations: 1.984 GiB, 3.84% gc time) # 1st ever without assigning
# 0.538783 seconds (104.39 k allocations: 494.558 MiB) # 2nd

# Soo CSV.read(file, DataFrame) seems way better only for 2nd loading of very
# big DataFrames