#Genotypic Similarity: linear kernel
using NamedArrays
function linear_GAMuT_geno(x)
    out = ["list", 2]
    names(out) = ["lc", "ev_Lc"]

    lc = x'*x #Made specifically for eigvals function
    ev_Lc = eigvals(lc, permute=false, scale=false)
    lc = A_mul_Bt(x, x)
    for ev in ev_Lc
        if ev > 1e-08
            ev = ev_Lc
        end
    end
    return([lc,ev_Lc])
end
