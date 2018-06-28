numtrios = 500
maf = 0.25

function minorAlleleCountsBinomial(numtrios, MAF)
    mom1bin = Binomial(1, MAF)
    mom1 = rand(mom1bin, numtrios)
    mom1[mom1 .== 0] = 2 
    mom2bin = Binomial(1, MAF)
    mom2 = rand(mom2bin, numtrios)
    mom2[mom2 .== 0] = 2
    
    dad1bin = Binomial(1, maf)
    dad1 = rand(dad1bin, numtrios)
    dad1[dad1 .== 0] = 2
    
    dad2bin = Binomial(1, maf)
    dad2 = rand(dad2bin, numtrios)
    dad2[dad2 .== 0] = 2
    
    eur = collect(hcat(mom1,mom2,dad1,dad2))
    

    kids1 = sample([1,2],numtrios, replace=true)
    kids2 = sample([3,4],numtrios, replace=true)


    kids = zeros(Int64, numtrios, 2)
    for id in 1:numtrios
        kids[id,1] = eur[id, kids1[id]]
        kids[id,2] = eur[id, kids2[id]] 
    end
    kids = collect(kids)
    eur_kid = (kids[:, 1] .== 1) + (kids[:, 2] .== 1)
    eur_mom = (eur[:, 1] .== 1) + (eur[:, 1] .== 1)
    eur_dad = (eur[:, 3] .== 1) + (eur[:, 4] .== 1)

    return(G_kid=eur_kid, G_mom=eur_mom, G_dad=eur_dad)
end