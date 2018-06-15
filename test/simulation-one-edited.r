library(MASS)
library(corpcor)

#variables
maf = 0.25 
n_unrelated = 5000 
n_variants = 500

minorAlleleCountsBinomial = function(numtrios,MAF){
    mom1 = rbinom(numtrios, 1, MAF)
    mom1[mom1==0]=2 ##to match cosi
    mom2 = rbinom(numtrios, 1, MAF)
    mom2[mom2==0]=2
    dad1 = rbinom(numtrios, 1, MAF)
    dad1[dad1==0]=2
    dad2 = rbinom(numtrios, 1, MAF)
    dad2[dad2==0]=2
    eur = cbind(mom1,mom2,dad1,dad2)

    kids1 <- sample(c(1,2), numtrios, replace=TRUE)
    kids2 <- sample(c(3,4), numtrios, replace=TRUE)

    kids <- matrix(numeric(numtrios*2), nrow=numtrios)
    for(id in 1:numtrios){
        kids[id,1] <- eur[id, kids1[id]]
        kids[id,2] <- eur[id, kids2[id]]
    }
    eur_kid <- (kids[,1] == 1) + (kids[,2] == 1)
    eur_mom <- (eur[,1] == 1) + (eur[,1] == 1)
    eur_dad <- (eur[,3] == 1) + (eur[,4] == 1)

    return( list(G_kid=eur_kid, G_mom=eur_mom, G_dad=eur_dad) )
}

##------------------------------------------------------------------------------
## Simulating genotypes:
##------------------------------------------------------------------------------

UNR_OBS = matrix(rep(NA,n_variants*n_unrelated),nrow=n_unrelated)
for(i in 1:n_variants){
    out = minorAlleleCountsBinomial(n_unrelated, maf)
    UNR_OBS[,i]= as.matrix(out$G_mom)
}