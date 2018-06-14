#GAMuT Test
using RCall
R"library(CompQuadForm)"
function testGAMuT(yc, λy, xc, λx )
    #Test statistic
    m = size(yc, 1)
    gamut = (1/m)*sum(sum(At_mul_B(yc, xc)))

    #Derive p-value of GAMuT statistic
    #Form vector of eigenvalue products 
    z = vec(A_mul_Bt([λy], [λx]))
    zsort = [sort(collect(Iterators.flatten(z)), rev=true)]
    scoredavies = [gamut*m^2]
    @rput zsort
    @rput scoredavies 
    zsortnum = R"as.numeric(unlist(zsort))"
    @rput zsortnum
    R"results_score = davies(scoredavies, zsortnum)"
    @rget results_score
    davies_pvalue = results_score[:Qq]
    return(davies_pvalue)
end
