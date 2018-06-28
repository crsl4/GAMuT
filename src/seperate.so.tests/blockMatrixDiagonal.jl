

function blockMatrixDiagonal(...)
    matrixList = [...]
    if matrixList[[1]] == [...]
        matrixList = matrixList[[1]]
        dimensions = [size(matrixList,1) ]
        finalDimension = sum(dimensions)
        finalMatrix = zeros(Int64, finalDimension, finalDimension)
        index = 1
        for k in 1:length(dimensions)
            finalMatrix[(index:(index+dimensions[k]-1),index:(index+dimensions[k]-1))] = matrixList[k]
            index=index+dimensions[k]
        end
        return(finalMatrix)
    end
end






