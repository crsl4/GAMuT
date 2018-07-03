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
        P0_UNR[i,:] = rand(MvNormal(mod(betamat_unr, G[i,causal_ind]), cov_unr),1)
    return(P0_UNR)
    end
end

simulatePhenotypes(2, "low", causal_ind, 1, "rare", MAF_unr, 100, rand(6, 6))