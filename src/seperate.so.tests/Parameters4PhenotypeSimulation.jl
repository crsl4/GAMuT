include("./src/seperate.so.tests/variables.jl")
include("./src/seperate.so.tests/createCovMatrix.jl")
using Rmath

## function to set up the parameters to simulate the phenotypes associated with genotypes
## it returns the beta matrix and the covariance matrix
function parameters4phenotypeSimulation(npheno, traitcor, causal_ind, nassoc, variant, MAF_unr)
    if npheno < nassoc
        error("Error: npheno<nassoc")
    end
    MAF_C_unr = MAF_unr[causal_ind]
    ## pairwise similarity randomly generated from unif(phencor_ll, phencor_ul)
    cov_unr= createCovMatrix(npheno, traitcor)
    ## beta matrix: npheno by num of causal variants
    betamat_unr = zeros(npheno, length(causal_ind))
    if nassoc > 0
        hvec_unr = fill(0.0, nassoc)
        for i in 1:nassoc
            if (variant=="rare")
                betamat_unr[i,:] = (0.4 + rnorm(length(causal_ind), 0, 0.1))*abs(log10.(MAF_C_unr))
            elseif (variant=="common")
                betamat_unr[i,:] = fill(log(1.5), length(causal_ind))
            end
            tmp = betamat_unr[i,:].^2
            tmp2 = 2*MAF_C_unr.*(1-MAF_C_unr)
            hvec_unr[i] = sum(tmp .* tmp2')
        end

            ## note: the first nassoc phenotypes are the ones that are associated with the genotype
        for i in 1:nassoc 
            for ii in 1:nassoc
                if i==ii
                    cov_unr[i,ii] = 1-hvec_unr[i]
                elseif ii>i 
                    cov_unr[i,ii] = cov_unr[i,ii]*((1-hvec_unr[i])^0.5)*((1-hvec_unr[ii])^0.5)
                    cov_unr[ii,i] = cov_unr[i,ii]
                end
            end
        end
    
        if !isposdef(cov_unr)
            isposdef!(cov_unr) || error("cannot make positive definite")
        end
    end
    return(betamat_unr, cov_unr)
end