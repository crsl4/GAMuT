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
            cor = rand(Uniform(phencor_ll,phencor_ul), Int(npheno*(npheno-1)/2))
            mat = zeros(npheno, npheno)
            #Lhttps://stackoverflow.com/questions/51068263/lower-triangular-matrix-equal-to-value-in-julia
            mat[tril!(trues(size(mat)), -1)] = cor
            #https://github.com/JuliaLang/julia/issues/23156
            mat = mat + mat' + Array(Diagonal(ones(npheno)))
        end
    end
    return(mat)
end
