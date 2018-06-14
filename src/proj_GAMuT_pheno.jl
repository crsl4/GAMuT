#Genotypic Similarity: Projection Kernel
using NamedArrays
using DataFrames
function proj_GAMuT_pheno(x)
    out = ["list", 2]
    names(out) = ["kc", "ev_Kc"]
    kc = x*inv(At_mul_B(x, x))*x'
    ev_Kc = repeat([1], inner = size(x, 2))
    out = [kc, ev_Kc]
    return(out)
end