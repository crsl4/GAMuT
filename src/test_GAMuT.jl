#GAMuT Test
R"library(CompQuadForm)"
function testGAMuT(yc, λy, xc, λx )
    #Test statistic
    m = size(yc, 1)
    gamut = (1/m)*sum(yc' .* xc) #Broadcast necessary for significant p-value
    
    #Derive p-value of GAMuT statistic
    #Form vector of eigenvalue products
    z = λy.*λx'
    #https://stackoverflow.com/questions/50822442/how-to-sort-a-matrix-in-julia
    zsort = [sort(collect(Iterators.flatten(z)), rev=true)]
    scoredavies = [gamut*m^2]
    @rput zsort
    @rput scoredavies
    R"zsortnum = as.numeric(unlist(zsort))"
    R"results_score = davies(scoredavies, zsortnum)"
    @rget results_score
    davies_pvalue = results_score[:Qq]
    return(davies_pvalue)
end
