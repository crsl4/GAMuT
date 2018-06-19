library(MASS)
library(corpcor)

##------------------------------------------------------------------------------
## Variables:  
##------------------------------------------------------------------------------
maf = 0.25 
traitcor = "med"
nassoc1 = 1 
nassoc2 = 1 
npheno1 = 1 
npheno2 = 1 
n_unrelated = 5000 
n_variants = 500
causal_var = 0.01
test.approach = 1 
sim.approach = 3 
ignoreZ = TRUE

##------------------------------------------------------------------------------
## Preliminary Functions: 
##------------------------------------------------------------------------------
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

## function to create the phenotype covariance matrix
## it takes as input the number of phenotypes: npheno, and the type of correlation:
## "none","low", "med", "high"
## "block" = 1-2, 3-4, "antiblock": only correlation btw phenotype matrices, not within
createCovMatrix = function(npheno, traitcor){
    if(traitcor == "block" || traitcor == "antiblock"){
        if(npheno %% 2 != 0)
            stop("cannot do block covariance matrix if npheno%%2 != 0")
    }
    if(traitcor == "none"){
        phencor_ll <- 0
        phencor_ul <- 0
    }else if (traitcor=="low") {
        phencor_ll <- 0
        phencor_ul <- 0.3
    } else {
        if (traitcor=="med" || traitcor == "block" || traitcor == "antiblock") {
            phencor_ll <- 0.3
            phencor_ul <- 0.5
        } else {
            if (traitcor=="high") {
                phencor_ll <- 0.5
                phencor_ul <- 0.7
            }
        }
    }
    if(npheno == 1){
        mat = as.matrix(runif(1,min=phencor_ll,max=phencor_ul))
    }else{
        if(traitcor != 'block' && traitcor != 'antiblock'){
            cor = runif((npheno*(npheno-1))/2,min=phencor_ll,max=phencor_ul)
            mat <- matrix(0, npheno,npheno)
            mat[lower.tri(mat, diag=FALSE)] <- cor
            mat = mat + t(mat) + diag(1, nrow=npheno)
        }else{
            if(traitcor == 'block'){
                m = npheno/2
                cor1 = runif((m*(m-1))/2,min=phencor_ll,max=phencor_ul)
                mat1 <- matrix(0, m,m)
                mat1[lower.tri(mat1, diag=FALSE)] <- cor1
                mat1 = mat1 + t(mat1) + diag(1, nrow=m)
                cor2 = runif((m*(m-1))/2,min=phencor_ll,max=phencor_ul)
                mat2 <- matrix(0, m,m)
                mat2[lower.tri(mat2, diag=FALSE)] <- cor2
                mat2 = mat2 + t(mat2) + diag(1, nrow=m)
                mat = blockMatrixDiagonal(mat1,mat2)
            }else{
                cor = runif((npheno*(npheno-1))/2,min=phencor_ll,max=phencor_ul)
                mat1 <- matrix(0, npheno,npheno)
                mat1[lower.tri(mat1, diag=FALSE)] <- cor
                mat1 = mat1 + t(mat1) + diag(1, nrow=npheno)
                mat = blockMatrixAntiDiagonal(mat1,npheno/2)
            }
        }
    }
    return(mat)
}

## function to set up the parameters to simulate the phenotypes associated with genotypes
## it returns the beta matrix and the covariance matrix
parameters4phenotypeSimulation = function(npheno, traitcor, causal.ind, nassoc, variant, MAF_unr){
    if(npheno < nassoc)
        stop("Error: npheno<nassoc")
    MAF_C_unr<-MAF_unr[causal.ind]
    ## pairwise similarity randomly generated from unif(phencor_ll, phencor_ul)
    cov_unr= createCovMatrix(npheno, traitcor)
    ## beta matrix: npheno by num of causal variants
    betamat_unr <- matrix(0,nrow=npheno,ncol=length(causal.ind))
    if(nassoc > 0){
        hvec_unr<-rep(0.0,nassoc)
        for (i in 1:nassoc) {
            if (variant=="rare")
                betamat_unr[i,] <- (0.4 + rnorm(length(causal.ind), 0, 0.1))*abs(log(MAF_C_unr, base=10))
            else if (variant=="common")
                betamat_unr[i,] <- rep(log(1.5),length(causal.ind))
            hvec_unr[i] <- sum(betamat_unr[i,]^2*2*MAF_C_unr*(1-MAF_C_unr))
        }

        ## note: the first nassoc phenotypes are the ones that are associated with the genotype
        for (i in 1:nassoc) {
            for (ii in 1:nassoc) {
                if (i==ii){
                    cov_unr[i,ii] <- 1-hvec_unr[i]
                } else if(ii>i) {
                    cov_unr[i,ii] <- cov_unr[i,ii]*((1-hvec_unr[i])^0.5)*((1-hvec_unr[ii])^0.5)
                    cov_unr[ii,i] <- cov_unr[i,ii]
                }
            }
        }

        if(!is.positive.definite(cov_unr))
            cov_unr <- make.positive.definite(cov_unr)
    }
    return( list(betamat=betamat_unr, cov=cov_unr) )
}

simulatePhenotypes = function(npheno, traitcor, causal.ind, nassoc, variant, MAF, n_unrelated, G){
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno, traitcor, causal.ind, nassoc, variant, MAF)
    betamat_unr = out2$betamat
    cov_unr = out2$cov
##    print(betamat_unr)
##    print(cov_unr)
##    print(nrow(UNR_OBS[1,causal.ind]))
##    print(ncol(UNR_OBS[1,causal.ind]))
    ## Actual phenotype simulation:
    P0_UNR <- matrix(numeric(n_unrelated*npheno), ncol=npheno)
    for(i in 1:n_unrelated)
        P0_UNR[i,] <- mvrnorm(1,betamat_unr %*% G[i,causal.ind], cov_unr)
    return( P0_UNR )
}

## Get Y1 and Y2 from phenotype matrix
splitPhenotypeMatrix = function(P,nassoc1,nassoc2, npheno1,npheno2){
    n = nrow(P)
    m = ncol(P)
    if(nassoc1>npheno1)
        stop("Error: nassoc1>npheno1")
    if(nassoc2>npheno2)
        stop("Error: nassoc2>npheno2")
    if(npheno1+npheno2 != m)
        stop("Error: npheno1+npheno2 != npheno")

    if(nassoc1 == 0){
        Y2 = P[,1:npheno2]
        Y1 = P[,(npheno2+1):m]
    }else if(nassoc2 == 0){
        Y1 = P[,1:npheno1]
        Y2 = P[,(npheno1+1):m]
    } else {
        Y1 = P[,c(1:nassoc1,(nassoc1+nassoc2+1):(nassoc1+nassoc2+(npheno1-nassoc1)))]
        Y2 = P[,c((nassoc1+1):(nassoc1+nassoc2),(nassoc1+nassoc2+npheno1-nassoc1+1):m)]
    }
    return(list(Y1=Y1,Y2=Y2))
}

## Each row in Y1 is simulated as normal with the row of Y2 as mean (assuming beta=1), and cov=I
simulateFullyMediatedPhenotype = function(Y2,npheno1,n_unrelated,traitcor){
    Y1 <- matrix(numeric(npheno1*n_unrelated), ncol=npheno1)
    cov = createCovMatrix(npheno1, traitcor)
    for(i in 1:n_unrelated){
        if( ncol(Y2) <= npheno1 ){
            mu = rep(0,npheno1)
            for(j in 1:ncol(Y2))
                mu[j] = Y2[i,j]
        }else {
            mu = Y2[i,1:npheno1]
        }
        Y1[i,] <- mvrnorm(1,mu, cov)
    }
    return( Y1 )
}

## Each row in Y1 is simulated as function of Y2 and G
simulatePartiallyMediatedPhenotype = function(Y2,npheno1,n_unrelated, traitcor,causal.ind,nassoc1,variant,MAF_unr,UNR_OBS){
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno1, traitcor, causal.ind, nassoc1, variant, MAF_unr)
    betamat_unr = out2$betamat
    cov_unr = out2$cov
    ## Actual phenotype simulation:
    Y1 <- matrix(numeric(n_unrelated*npheno1), ncol=npheno1)
    for(i in 1:n_unrelated){
        if( ncol(Y2) <= npheno1 ){
            mu = rep(0,npheno1)
            for(j in 1:ncol(Y2))
                mu[j] = Y2[i,j]
        }else {
            mu = Y2[i,1:npheno1]
        }
        Y1[i,] <- mvrnorm(1,betamat_unr %*% UNR_OBS[i,causal.ind]+mu, cov_unr)
    }
    return( Y1 )
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
