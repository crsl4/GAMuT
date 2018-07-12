#Genotypic Similarity: linear kernel
function linear_GAMuT_geno(x)
    lc = x'*x #Made specifically for eigvals function
    ev_Lc = eigvals(lc, permute=false, scale=false)
    lc = A_mul_Bt(x, x)
    ev_Lc = ev_Lc[ev_Lc .> 1e-08]
    sort!(ev_Lc,rev=true)
    return([lc,ev_Lc])
end
