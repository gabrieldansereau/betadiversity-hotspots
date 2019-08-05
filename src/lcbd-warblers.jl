using Distributed
addprocs(9)

@time @everywhere include("src/required.jl")

## Get & prepare data
@time @everywhere begin
    # Load data from CSV files
    df = CSV.read("../data/warblers_qc_2018.csv", header=true, delim="\t")
    # Prepare data (select columns, arrange values)
    df = prepare_csvdata(df)
    # Separate species
    taxa_occ = [df[df.species .== u,:] for u in unique(df.species)]

    # Define coordinates range
    lon_range = (-136.0, -58.0)
    lat_range = (40.5, 56.0)
end

## Get the worldclim data
@time wc_vars = pmap(x -> worldclim(x)[lon_range, lat_range], 1:19);

# Make predictions for all species
@time predictions = pmap(x -> species_bclim(x, wc_vars), taxa_occ);

## Get the LCBD
begin
    # Create Y -> site-by-species community data table
    Y = zeros(Int64, (prod(size(predictions[1])),length(taxa_occ)))
    # Fill Y with community predictions
    @progress for gc in eachindex(predictions[1].grid) # loop for all sites
        # Group predictions for all species in site
        R = map(x -> x.grid[gc], predictions)
        # Fill Y with binary values -> 1 if species prediction for site != NaN, 0 if == NaN
        global Y[gc,:] = .!isnan.(R)
    end

    ## Compute statistics
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

    ## Arrange LCBD values as grid
    # Create empty grid
    t_lcbd = zeros(Float64, size(predictions[1]))
    # Scale LCBDi values to maximum value
    LCBDi = LCBDi./maximum(LCBDi)
    # Fill-in grid
    for i in eachindex(t_lcbd)
        # Add LCBD values only if prediction != NaN
        t_lcbd[i] = any(Y[i,:] .> 0) ? LCBDi[i] : NaN # ?: is ternary operator, ~if-else
    end
    #=
    ## Show error in initial code
    # Reverse column order (last colimn/species has very few obs)
    Y = Y[:,sort(1:end, rev=true)]
    # Fill-in grid
    for i in eachindex(t_lcbd)
        # Add LCBD values only if prediction != NaN
        t_lcbd[i] = Y[i] > 0 ? LCBDi[i] : NaN # ?: is ternary operator, ~if-else
    end
    =#

    # Create SDMLayer with LCBD values
    LCBD = SDMLayer(t_lcbd, predictions[1].left, predictions[1].right, predictions[1].bottom, predictions[1].top)
end

## Plot results
sdm_plot = plotSDM(LCBD, type="lcbd")

## Save result
savefig(sdm_plot, "fig/warblers/lcbd-map-qc2018.pdf")
