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
