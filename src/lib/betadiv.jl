#### Beta diversity calculation functions

## Function to calculate beta diversity statistics
function betadiv(Y)
    # S -> squared deviations from column mean
    S = (Y .- mean(Y; dims=1)) .^ 2.0
    # SStotal -> total sum of squares
    SStotal = sum(S)
    # BDtotal -> index of beta diversity, unbiased & comparable estimator of Var(Y)
    BDtotal = SStotal / (size(Y, 1) - 1)
    # SSj -> sum of squares for species j
    SSj = sum(S; dims=1)
    # SCBDj -> species contribution to beta diversity (species j, relative)
    SCBDj = SSj ./ SStotal
    # SSi -> sum of squares for site i
    SSi = sum(S; dims=2)
    # LCBD -> local contribution to beta diversity (site i, relative)
    LCBDi = SSi ./ SStotal
    # Combine results in tuple
    res = (
        S=S, SStotal=SStotal, BDtotal=BDtotal, SSj=SSj, SCBDj=SCBDj, SSi=SSi, LCBDi=LCBDi
    )
    return res
end

## Function for permutation tests
function permtest(Y, res)
    # Permutation of matrix Y
    Y_perm = hcat(shuffle.(eachcol(Y))...)
    # Recalculate BD statistics
    res_p = betadiv(Y_perm)
    # Test if permuted LCBD is greater than original LCBD
    ge = res_p.LCBDi .>= res.LCBDi
    return ge
end

## Function combining BD calculation & permutation tests
function betadiv_perm(Y; nperm=999, distributed=true)
    n = size(Y, 1)
    p = size(Y, 2)
    # Initial BD results
    res = betadiv(Y)
    # Permutations
    if nperm > 0
        nGE_L = ones(Int64, n)
        # Permutation test, clumsy parallelization
        ge = @showprogress pmap(x -> permtest(Y, res), 1:nperm; distributed=distributed)
        # Compile number of permuted LCBDs greater than original LCBD
        for bitarray in ge
            nGE_L[findall(bitarray)] .+= 1
        end
        # Standardize counts
        p_LCBD = nGE_L ./ (nperm + 1)
    else
        p_LCBD = nothing
    end
    # Combine results in tuple
    res_perm = (res..., pLCBD=p_LCBD)
    return res_perm
end
