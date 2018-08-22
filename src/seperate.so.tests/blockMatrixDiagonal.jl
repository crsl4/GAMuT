## Does not work in 0.6 or 0.7 (ellipsis) 
## This is not used in simulation and may be for a specific case scenario. 

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






