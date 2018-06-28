using Distributions
##------------------------------------------------------------------------------
## Variables:  
##------------------------------------------------------------------------------
maf = 0.25 
traitcor = "med"
nassoc1 = 1 
nassoc2 = 1 
npheno1 = 1 
npheno2 = 1 
n_unrelated = 5000 
n_variants = 500
causal_var = 0.01
test_approach = 1 
sim_approach = 3 
ignoreZ = true

function simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF, n_unrelated, G)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno, traitcor, causal_ind, nassoc, variant, MAF)
    betamat_unr = out2[betamat]
    cov_unr = out2[cov]
    ## Actual phenotype simulation:
    P0_UNR = [convert(Int64, n_unrelated*npheno), npheno]
    for i in 1:n_unrelated
        P0_UNR[i,] = MvNormal(mod(betamat_unr, G[i,causal_ind]), cov_unr)
    return(P0_UNR)
    end
end