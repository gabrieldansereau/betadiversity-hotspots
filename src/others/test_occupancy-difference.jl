using Pkg: Pkg
Pkg.activate(".")
using Test
include("../required.jl")

# Checkout 41a9e9 & load data
verify_jld2_data("data/jld2"; extract_recent=true)
@load joinpath("data", "jld2", "raw-distributions.jld2") distributions spenames specommon speindex
Y = calculate_Y(distributions)
inds_obs = _indsobs(Y)
Yobs = _Yobs(Y, inds_obs)
env_df = CSV.read(
    joinpath("data", "proc", "distributions_env_full.csv"),
    DataFrame;
    header=true,
    delim="\t",
)
# Save new data
new = (
    distributions=distributions,
    spenames=spenames,
    Y=Y,
    inds_obs=inds_obs,
    Yobs=Yobs,
    env_df=env_df,
)
# Checkout a1cbac, re-run, & save old data
old = (
    distributions=distributions,
    spenames=spenames,
    Y=Y,
    inds_obs=inds_obs,
    Yobs=Yobs,
    env_df=env_df,
)

# Investigate sites with NA for land cover
allowmissing!(old.env_df)
for col in eachcol(old.env_df)
    replace!(col, NaN => missing)
end
@test nrow(old.env_df[old.inds_obs, :]) == nrow(dropmissing(old.env_df[old.inds_obs, :]))

# Investigate number of occurrence differece
old.distributions
new.distributions
comp = DataFrame(; old=length.(old.distributions), new=length.(new.distributions))
insertcols!(comp, :diff => comp.new .- comp.old)
insertcols!(comp, :reldiff => comp.diff ./ comp.old)
show(comp; allrows=true)
# Verify order of species
@test old.spenames == new.spenames
insertcols!(comp, :oldnames => old.spenames, :newnames => new.spenames)
insertcols!(comp, :equal => isequal.(comp.oldnames, comp.newnames))
show(comp; allrows=true)

# Compare by species
olddf = DataFrame(; sp=old.spenames, n_old=length.(old.distributions))
newdf = DataFrame(; sp=new.spenames, n_new=length.(new.distributions))
comp2 = outerjoin(olddf, newdf; on=:sp)
insertcols!(comp2, :diff => comp2.n_new .- comp2.n_old)
insertcols!(comp2, :reldiff => comp2.diff ./ comp2.n_old)
show(comp2; allrows=true)
show(sort(comp2, :diff); allrows=true)
show(sort(comp2, :reldiff); allrows=true)
