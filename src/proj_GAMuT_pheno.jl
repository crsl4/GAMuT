#Genotypic Similarity: Projection Kernel
function proj_GAMuT_pheno(x)
    kc = x*inv(At_mul_B(x, x))*x'
    ev_Kc = repeat([1], inner = size(x, 2))
    out = [kc, ev_Kc]
    return(out)
end
