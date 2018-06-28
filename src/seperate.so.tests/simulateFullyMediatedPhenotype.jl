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

function simulateFullyMediatedPhenotype(Y2,npheno1,n_unrelated,traitcor)
    Y1 = zeros(npheno1*n_unrelated, npheno1)
    cov = createCovMatrix(npheno1, traitcor)
    for i in 1:n_unrelated
        if size(Y2, 2) <= npheno1
            mu = repeat(0,outer=npheno1)
            for j in 1:size(Y2, 2)
                mu[j] = Y2[i,j]
            end
        else 
            mu = Y2[i,1:npheno1]
        end
        d = MvNormal(mu, cov)
        Y1[i,] = rand(d,1)
    end
    return(Y1)
end