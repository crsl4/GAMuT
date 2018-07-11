using Distributions

include("./src/seperate.so.tests/variables.jl")
include("./src/seperate.so.tests/createCovMatrix.jl")

function simulateFullyMediatedPhenotype(Y2,npheno1,n_unrelated,traitcor)
    Y1 = zeros(npheno1*n_unrelated, npheno1)
    cov = createCovMatrix(npheno1, traitcor)
    for i in 1:n_unrelated
        if size(Y2, 2) <= npheno1
            mu = zeros(Int64, npheno1)
            for j in 1:size(Y2, 2)
                mu[j] = Y2[i,j]
            end
        else  #Works  
            mu = Y2[i,1:npheno1]
        end
        d = MvNormal(mu, cov)
        Y1[i,:] = rand(d,1)
    end
    return(Y1)
end

npheno1=2
simulateFullyMediatedPhenotype(Y2,npheno1,n_unrelated,traitcor)
###InexactError