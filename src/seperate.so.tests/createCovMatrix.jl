using Distributions

npheno = 1
traitcor = "med"

## function to create the phenotype covariance matrix
## it takes as input the number of phenotypes: npheno, and the type of correlation:
## "none","low", "med", "high"
## "block" = 1-2, 3-4, "antiblock": only correlation btw phenotype matrices, not within
function createCovMatrix(npheno, traitcor)
    if traitcor == "block" || traitcor == "antiblock"
        if(mod(npheno, 2) != 0)
            error("Error: Cannot do block covariance matrix if npheno%%2 != 0")
        end
    end
    if traitcor == "none"
        phencor_ll = 0
        phencor_ul = 1e-12
    elseif traitcor == "low"
        phencor_ll = 0
        phencor_ul = 0.3
    else
        if traitcor == "med" || traitcor == "block" || traitcor == "antiblock"
            phencor_ll = 0.3
            phencor_ul = 0.5
         else 
            if traitcor == "high"
                phencor_ll = 0.5
                phencor_ul = 0.7
            end
        end
    end

    if npheno == 1
        mat = rand(Uniform(phencor_ll,phencor_ul),1,1)
    else
        if traitcor != "block" && traitcor != "antiblock"
            cor = rand(Uniform(phencor_ll,phencor_ul), size(npheno*(npheno-1)/2,1), 1)
            mat = zeros(npheno, npheno)
            #Lhttps://stackoverflow.com/questions/51068263/lower-triangular-matrix-equal-to-value-in-julia
            mat[tril!(trues(size(mat)), -1)] = cor
            mat = mat + mat' + eye(npheno)
        else
            if traitcor == "block"
                m = npheno/2
                cor1 = rand(Uniform(phencor_ll,phencor_ul), size(m*(m-1)/2,1), 1)
                mat1 = zeros(m, m)
                cor1 = mat1[tril!(trues(size(mat1)), -1)] = cor1
                mat1 = mat1 + mat1' + eye(m)
                cor2 = rand(Uniform(phencor_ll,phencor_ul), size(m*(m-1)/2,1), 1)
                mat2 = [0, m,m]
                cor2 = mat2[tril!(trues(size(mat2)), -1)] = cor2
                mat2 = mat2 + mat2' + eye(m)
                mat = blockMatrixDiagonal(mat1, mat2)
            else
                cor = rand(Uniform(phencor_ll,phencor_ul), size(npheno*(npheno-1)/2), 1)
                mat1 = [0, npheno, npheno]
                cor = mat1[tril!(trues(size(mat1)), -1)] = cor
                mat1 = mat1 + mat1' + eye(npheno)
                mat = blockMatrixAntiDiagonal(mat1,npheno/2)
            end
        end
    end
    return(mat)
end
