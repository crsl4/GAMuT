library(MASS)
library(corpcor)

#variables
maf = 0.25 
traitcor = med 
nassoc1 = 1 
nassoc2 = 1 
npheno1 = 1 
npheno2 = 1 
n_unrelated = 5000 
n_variants = 500
test_approach = 1 
simulation_approach = 3 
ignore_Z = TRUE

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

simulatePhenotypesMediation = function(npheno1, npheno2,traitcor, nassoc1, nassoc2, causal.ind, MAF_unr, n_unrelated, variant,UNR_OBS, approach=1, numpcs=1, ignoreZ=FALSE){
    if(approach == 1){
        P0_UNR = simulatePhenotypes(npheno1+npheno2, traitcor, causal.ind, nassoc1+nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
        ## Now, we want to split this matrix into two to test for mediation
        out3 = splitPhenotypeMatrix(P0_UNR,nassoc1, nassoc2,npheno1,npheno2)
        Y1 = out3$Y1
        Y2 = out3$Y2
    }else if(approach == 20){ ##approach2: no mediation
        Y1 = simulatePhenotypes(npheno1, traitcor, causal.ind, nassoc1, variant, MAF_unr, n_unrelated, UNR_OBS)
        Y2 = simulatePhenotypes(npheno2, traitcor, causal.ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
    }else if(approach == 2 && nassoc1 == 0){ ##approach2: full mediation
        Y2 = simulatePhenotypes(npheno2, traitcor, causal.ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
        Y1 = simulateFullyMediatedPhenotype(Y2, npheno1, n_unrelated, traitcor)
    }else if(approach == 2 && nassoc1 > 0){ ##approach2: partial mediation
        Y2 = simulatePhenotypes(npheno2, traitcor, causal.ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
        Y1 = simulatePartiallyMediatedPhenotype(Y2,npheno1,n_unrelated, traitcor,causal.ind,nassoc1,variant,MAF_unr,UNR_OBS)
    }else if(approach == 3){
        Z = as.matrix(rnorm(n_unrelated))
        if(nassoc2 == 0){ ##approach3: no association with G
            if(nassoc1 != 0)
                stop("nassoc1 cannot be != 0 if nassoc2 is ==0")
            P2 = simulateFullyMediatedPhenotype(Z, npheno2, n_unrelated, traitcor)
            P1 = simulateFullyMediatedPhenotype(Z, npheno1, n_unrelated, traitcor)
        }else if(nassoc2 > 0){ ##approach3: Y2 associated with G
            P2 = simulatePartiallyMediatedPhenotype(Z,npheno2,n_unrelated, traitcor,causal.ind,nassoc2,variant,MAF_unr,UNR_OBS)
            if(nassoc1 == 0){
                P1 = simulateFullyMediatedPhenotype(Z, npheno1, n_unrelated, traitcor)
            }else if(nassoc1 > 0){
                P1 = simulatePartiallyMediatedPhenotype(Z,npheno1,n_unrelated, traitcor,causal.ind,nassoc1,variant,MAF_unr,UNR_OBS)
            }else{
                stop("nassoc1 should be >=0")
            }
        }
        if(ignoreZ){
            Y1 = P1
            Y2 = P2
        }else{
            ## PCA ---------------------------------------
            P = cbind(P1,P2)
            pc = prcomp(P, scale. = TRUE, center = TRUE)
            prop = pc$sdev[1]/sum(pc$sdev)
            print(paste("Proportion of variance explained by PC1:",prop))
            Z = as.matrix(pc$x[,1:numpcs])
            ## -------------------------------------------
            Y2 = matrix(numeric(nrow(P2)*ncol(P2)), ncol=ncol(P2))
            for(i in 1:ncol(P2)){
                f = lm(P2[,i]~Z)
                Y2[,i] = residuals(f)
            }
            Y1 = matrix(numeric(nrow(P1)*ncol(P1)), ncol=ncol(P1))
            for(i in 1:ncol(P1)){
                f = lm(P1[,i]~Z)
                Y1[,i] = residuals(f)
            }
        }
    }
    return( list(Y1=Y1, Y2=Y2))
}


##------------------------------------------------------------------------------
## Simulating genotypes:
##------------------------------------------------------------------------------

UNR_OBS = matrix(rep(NA,n_variants*n_unrelated),nrow=n_unrelated)
for(i in 1:n_variants){
    out = minorAlleleCountsBinomial(n_unrelated, maf)
    UNR_OBS[,i]= as.matrix(out$G_mom)
}

n_causal = floor(causal_var*n_variants)
causal.ind = sample(1:n_variants,n_causal,replace=FALSE)

## determine the MAF of each variant in the sample
MAF_unr = colMeans(UNR_OBS)/2

## we set "variant" to use the same functions as with cosi
if(maf > 0.05){
    variant = 'common'
}else{
    variant = 'rare'
}


##------------------------------------------------------------------------------
## Simulating phenotypes:
##------------------------------------------------------------------------------
out2 = simulatePhenotypesMediation(npheno1, npheno2,traitcor, nassoc1, nassoc2, causal.ind, MAF_unr, n_unrelated, variant, UNR_OBS,approach=sim.approach, ignoreZ=ignoreZ)
Y1 = out2$Y1
Y2 = out2$Y2
