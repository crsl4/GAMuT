## Get Y1 and Y2 from phenotype matrix
function splitPhenotypeMatrix(P,nassoc1,nassoc2, npheno1,npheno2)
    n = size(P, 1)
    m = size(P, 2)
    if nassoc1>npheno1
        error("Error: nassoc1>npheno1")
    end
    if nassoc2>npheno2
        error("Error: nassoc2>npheno2")
    end
    if npheno1+npheno2 != m
        error("Error: npheno1+npheno2 != npheno")
    end
    if nassoc1 == 0
        Y2 = P[1:size(P,1), 1:npheno2]
        Y1 = P[1:size(P,1), (npheno2+1):m]
    elseif nassoc2 == 0
        Y1 = P[1:size(P, 1), 1:npheno1]
        Y2 = P[1:size(P, 1),(npheno1+1):m]
    else
        Y1 = P[:, vcat(1:nassoc1,(nassoc1+nassoc2+1):(nassoc1+nassoc2+(npheno1-nassoc1)))]
        Y2 = P[:, vcat((nassoc1+1):(nassoc1+nassoc2),(nassoc1+nassoc2+npheno1-nassoc1+1):m)]
    end
    return([Y1,Y2])
end

