using Distributions
using Rmath
using Distributions
using MultivariateStats
using GLM
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

function simulatePartiallyMediatedPhenotype(Y2,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno1, traitcor, causal_ind, nassoc1, variant, MAF_unr)
    betamat_unr = out2[:betamat]
    cov_unr = out2[:cov]
    ## Actual phenotype simulation:
    Y1 = zeros(npheno1*n_unrelated, npheno1)
    for i in 1:n_unrelated
        if size(Y2, 2) <= npheno1
            mu = repeat(0, outer=npheno1)
            for j in 1:size(Y2, 2)
                mu[j] = Y2[i, j]
            end
        else
            mu = Y2[i, 1:npheno1]
        end
        Y1[i,] = MvNormal(mod(betamat_unr, UNR_OBS[i,causal_ind]+mu), cov_unr)
    end
    return(Y1)
end