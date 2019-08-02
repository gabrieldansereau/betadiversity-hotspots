using CSV
using Random
using DataFrames
using Statistics

Y = CSV.read("../data/Y-can.csv", delim="\t")
Y = Matrix(Y)

function BD(Y)
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

res = BD(Y)

function BDperm(Y; nperm=999)
    n = size(Y, 1)
    p = size(Y, 2)
    # Initial BD results
    res = BD(Y)
    # Permutations
    if nperm > 0
        nGE_L = ones(Int64, n)
        for iperm in 1:nperm
            Y_perm = hcat(shuffle.(eachcol(Y))...)
            res_p = BD(Y_perm)
            ge = res_p.LCBDi .>= res.LCBDi
            nGE_L[findall(ge)] .+= 1
        end
        p_LCBD = nGE_L./(nperm+1)
    else
        p_LCBD = NaN
    end
    res_perm = (res..., pLCBD = p_LCBD)
    return res_perm
end
res_perm = BDperm(Y, nperm=0)
res_perm = BDperm(Y, nperm=99)

findall(res_perm.pLCBD .<= 0.05)
