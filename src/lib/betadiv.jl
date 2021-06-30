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
