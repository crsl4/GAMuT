#Genotypic Similarity: linear kernel
function linear_GAMuT_geno(x)
    lc = x'*x #Made specifically for eigvals function
    @rput lc
    R"""
    ev_Lc = eigen(lc, symmetric=TRUE, only.values=TRUE)$values
    """
    @rget ev_Lc
    lc = A_mul_Bt(x, x)
    ev_Lc = ev_Lc[ev_Lc .> 1e-08]
    sort!(ev_Lc,rev=true)
    return([lc,ev_Lc])
end
