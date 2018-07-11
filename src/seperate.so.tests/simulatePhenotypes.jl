using Distributions

include("./src/seperate.so.tests/variables.jl")
include("./src/seperate.so.tests/Parameters4PhenotypeSimulation.jl")

function simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF_unr, n_unrelated, G)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno, traitcor, causal_ind, nassoc, variant, MAF_unr)
    betamat_unr = out2[1]
    cov_unr = out2[2]
    ## Actual phenotype simulation:
    P0_UNR = [convert(Int64, n_unrelated*npheno), npheno]
    for i in 1:n_unrelated
        for x in size(causal_ind,2)
        P0_UNR[i,:] = rand(MvNormal(mod(betamat_unr, G[i,causal_ind[x]]), cov_unr),1)
        end
    return(P0_UNR)
    end
end

simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF_unr, n_unrelated, rand(100,2))






##What is G? currently is a 100 X 2
##ERROR: BoundsError: attempt to access 100×2 Array{Float64,2} at index [1, 261]