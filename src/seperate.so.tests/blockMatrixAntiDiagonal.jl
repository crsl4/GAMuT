## Works correctly in 0.6 but have not yet checking in 0.7
## This is not used in simulation and may be for a specific case scenario. 

## input only one matrix (and size of block) that will have the block diagonal
## replaced by identities
function blockMatrixAntiDiagonal(M,n)
    ntotal = size(M,1)
    if mod(ntotal, n) != 0
        error("Error: dimension of matrix M should be a multiple of n")
    end
    numI = ntotal / n
    #https://stackoverflow.com/questions/26042691/declaring-multiple-arrays-in-julia
    matrixList = [Array{Int64}(n, n) for i = 1:numI]
    for i in range(1, convert(Int64, numI))
        matrixList[i] = eye(n)
    end

dimensions = [size(matrixList[i],1) for i = 1:convert(Int64, numI)]
index = 1
    for k in range(1, length(dimensions))
        M[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)] = matrixList[k]
        index = index + dimensions[k]
    end
    return(M)
end