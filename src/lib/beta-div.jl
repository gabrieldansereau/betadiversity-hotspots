#### Beta diversity calculation functions
@everywhere using ProgressMeter

## Function to calculate beta diversity statistics
@everywhere function BD(Y)
    # S -> squared deviations from column mean
    S = (Y .- mean(Y; dims=1)).^2.0
    # SStotal -> total sum of squares
    SStotal = sum(S)
    # BDtotal -> index of beta diversity, unbiased & comparable estimator of Var(Y)
    BDtotal = SStotal / (size(Y,1)-1)
    # SSj -> sum of squares for species j
    SSj = sum(S; dims=1)
    # SCBDj -> species contribution to beta diversity (species j, relative)
    SCBDj = SSj ./ SStotal
    # SSi -> sum of squares for site i
    SSi = sum(S; dims=2)
    # LCBD -> local contribution to beta diversity (site i, relative)
    LCBDi = SSi ./ SStotal
    # Combine results in tuple
    res = (S = S, SStotal = SStotal, BDtotal = BDtotal,
            SSj = SSj, SCBDj = SCBDj, SSi = SSi, LCBDi = LCBDi)
    return res
end

## Function for permutation tests
@everywhere function permtest(Y, res)
    # Permutation of matrix Y
    Y_perm = hcat(shuffle.(eachcol(Y))...)
    # Recalculate BD statistics
    res_p = BD(Y_perm)
    # Test if permuted LCBD is greater than original LCBD
    ge = res_p.LCBDi .>= res.LCBDi
    return ge
end

## Function combining BD calculation & permutation tests
@everywhere function BDperm(Y; nperm=999, distributed=true)
    n = size(Y, 1)
    p = size(Y, 2)
    # Initial BD results
    res = BD(Y)
    # Permutations
    if nperm > 0
        nGE_L = ones(Int64, n)
        # Permutation test, clumsy parallelization
        ge = @showprogress pmap(x -> permtest(Y, res), 1:nperm, distributed=distributed)
        # Compile number of permuted LCBDs greater than original LCBD
        for bitarray in ge
            nGE_L[findall(bitarray)] .+= 1
        end
        # Standardize counts
        p_LCBD = nGE_L./(nperm+1)
    else
        p_LCBD = NaN
    end
    # Combine results in tuple
    res_perm = (res..., pLCBD = p_LCBD)
    return res_perm
end

########

## Test functions

#=
using Distributed
addprocs(9)
@everywhere begin
    using CSV
    using Random
    using DataFrames
    using Statistics
end

# Input Y matrix
Y = CSV.read("../data/Y-can.csv", delim="\t")
Y = Matrix(Y)

@time res = BD(Y)
@time permtest(Y, res)
@time res_perm = BDperm(Y, nperm=0);
@time res_perm = BDperm(Y, nperm=49);
# 5 sec distributed vs 15 sec not distributed
@time res_perm = BDperm(Y, nperm=999, distributed=false);
# 90 sec distributed vs 300 sec
# Parallelized function seems ~3x faster

# Find sites with significant LCBD contributions
findall(res_perm.pLCBD .<= 0.05)
=#
