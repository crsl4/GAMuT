using Distributions
using Rmath
using Distributions
using MultivariateStats
using GLM

include("./src/seperate.so.tests/variables.jl")
include("./src/seperate.so.tests/parameters4phenotypeSimulation.jl")

function simulatePartiallyMediatedPhenotype(Y2,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno1, traitcor, causal_ind, nassoc1, variant, MAF_unr)
    betamat_unr = out2[1]
    cov_unr = out2[2]
    ## Actual phenotype simulation:
    Y1 = zeros(n_unrelated, npheno1)
    for i in 1:n_unrelated
        if size(Y2, 2) <= npheno1
            mu = zeros(npheno1)
            for j in 1:size(Y2,2)
                mu = Y2[i, j]
            end  
        else
            #Works
            mu = Y2[i, 1:npheno1]
        end
        mu0 = betamat_unr * UNR_OBS[i,causal_ind]'
        mu0 = vcat(mu0...)
        Y1[i,:] = rand(MvNormal( mu0 + mu, cov_unr),1)
    end
    return(Y1)
end

simulatePartiallyMediatedPhenotype(Y2, 3, n_unrelated, traitcor, causal_ind, nassoc1, variant, MAF_unr, UNR_OBS)
##ERROR: DimensionMismatch("tried to assign 1 elements to 3 destinations") line 21

simulatePartiallyMediatedPhenotype(Y2, 2, n_unrelated, traitcor, causal_ind, nassoc1, variant, MAF_unr, UNR_OBS)
##ERROR: MethodError: no method matching Distributions.MvNormal(::Array{Float64,2}, ::Array{Float64,2})