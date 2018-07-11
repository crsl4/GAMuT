using Distributions
using MultivariateStats
using GLM
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
test_approach = 1 
sim_approach = 3 
ignoreZ = true
srand(1234)

##------------------------------------------------------------------------------
## Preliminary Functions: 
##------------------------------------------------------------------------------
function minorAlleleCountsBinomial(numtrios, maf)
    mom1bin = Binomial(1, maf)
    mom1 = rand(mom1bin, numtrios)
    mom1[mom1 .== 0] = 2 
    mom2bin = Binomial(1, maf)
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

    return(eur_kid, eur_mom, eur_dad)
end

# builds a block matrix whose diagonals are the square matrices provided.
# m1=matrix(runif(10*10),nrow=10,ncol=10)
# m2=matrix(runif(5*5),nrow=5,ncol=5)
# blockMatrix<-blockMatrixDiagonal(m1,m2,m2,m1)
# or
# blockMatrix<-blockMatrixDiagonal(list(m1,m2,m2,m1))
# C.Ladroue
# http://chrisladroue.com/2011/04/block-diagonal-matrices-in-r/
function blockMatrixDiagonal(...)
    matrixList = [...]
    if matrixList[[1]] == [...]
        matrixList = matrixList[[1]]
        dimensions = [size(matrixList,1) ]
        finalDimension = sum(dimensions)
        finalMatrix = zeros(Int64, finalDimension, finalDimension)
        index = 1
        for k in 1:length(dimensions)
            finalMatrix[(index:(index+dimensions[k]-1),index:(index+dimensions[k]-1))] = matrixList[k]
            index=index+dimensions[k]
        end
        return(finalMatrix)
    end
end

## input only one matrix (and size of block) that will have the block diagonal
## replaced by identities
function blockMatrixAntiDiagonal(M,n)
    ntotal = size(M,1)
    if mod(ntotal, n) != 0
        error("Error: dimension of matrix M should be a multiple of n")
    end
    numI = ntotal / n
    #https://stackoverflow.com/questions/26042691/declaring-multiple-arrays-in-julia
    matrixList = [Array{Int64}(n, n) for i = 1:numI]
    for i in range(1, convert(Int64, numI))
        matrixList[i] = eye(n)
    end

dimensions = [size(matrixList[i],1) for i = 1:convert(Int64, numI)]
index = 1
    for k in range(1, length(dimensions))
        M[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)] = matrixList[k]
        index = index + dimensions[k]
    end
    return(M)
end



## function to create the phenotype covariance matrix
## it takes as input the number of phenotypes: npheno, and the type of correlation:
## "none","low", "med", "high"
## "block" = 1-2, 3-4, "antiblock": only correlation btw phenotype matrices, not within
function createCovMatrix(npheno, traitcor)
    if traitcor == "block" || traitcor == "antiblock"
        if(mod(npheno, 2) != 0)
            error("Error: Cannot do block covariance matrix if npheno%%2 != 0")
        end
    end
    if traitcor == "none"
        phencor_ll = 0
        phencor_ul = 1e-12
    elseif traitcor == "low"
        phencor_ll = 0
        phencor_ul = 0.3
    else
        if traitcor == "med" || traitcor == "block" || traitcor == "antiblock"
            phencor_ll = 0.3
            phencor_ul = 0.5
         else 
            if traitcor == "high"
                phencor_ll = 0.5
                phencor_ul = 0.7
            end
        end
    end

    if npheno == 1
        mat = rand(Uniform(phencor_ll,phencor_ul),1,1)
    else
        if traitcor != "block" && traitcor != "antiblock"
            cor = rand(Uniform(phencor_ll,phencor_ul), size(npheno*(npheno-1)/2,1), 1)
            mat = zeros(npheno, npheno)
            #Lhttps://stackoverflow.com/questions/51068263/lower-triangular-matrix-equal-to-value-in-julia
            mat[tril!(trues(size(mat)), -1)] = cor
            mat = mat + mat' + eye(npheno)
        end
    end
    return(mat)
end



## function to set up the parameters to simulate the phenotypes associated with genotypes
## it returns the beta matrix and the covariance matrix
using Rmath
function parameters4phenotypeSimulation(npheno, traitcor, causal_ind, nassoc, variant, MAF_unr)
    if npheno < nassoc
        error("Error: npheno<nassoc")
    end
    MAF_C_unr = MAF_unr[causal_ind]
    ## pairwise similarity randomly generated from unif(phencor_ll, phencor_ul)
    cov_unr= createCovMatrix(npheno, traitcor)
    ## beta matrix: npheno by num of causal variants
    betamat_unr = zeros(npheno, length(causal_ind))
    if nassoc > 0
        hvec_unr = fill(0.0, nassoc)
        for i in 1:nassoc
            if (variant=="rare")
                betamat_unr[i,:] = (0.4 + rnorm(length(causal_ind), 0, 0.1))*abs(log10.(MAF_C_unr))
            elseif (variant=="common")
                betamat_unr[i,:] = fill(log(1.5), length(causal_ind))
            end
            tmp = betamat_unr[i,:].^2
            tmp2 = 2*MAF_C_unr.*(1-MAF_C_unr)
            hvec_unr[i] = sum(tmp .* tmp2')
        end

            ## note: the first nassoc phenotypes are the ones that are associated with the genotype
        for i in 1:nassoc 
            for ii in 1:nassoc
                if i==ii
                    cov_unr[i,ii] = 1-hvec_unr[i]
                elseif ii>i 
                    cov_unr[i,ii] = cov_unr[i,ii]*((1-hvec_unr[i])^0.5)*((1-hvec_unr[ii])^0.5)
                    cov_unr[ii,i] = cov_unr[i,ii]
                end
            end
        end
    
        if !isposdef(cov_unr)
            isposdef!(cov_unr) || error("cannot make positive definite")
        end
    end
    return(betamat_unr, cov_unr)
end


function simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF, n_unrelated, G)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno, traitcor, causal_ind, nassoc, variant, MAF)
    betamat_unr = out2[betamat]
    cov_unr = out2[cov]
    ## Actual phenotype simulation:
    P0_UNR = [convert(Int64, n_unrelated*npheno), npheno]
    for i in 1:n_unrelated
        P0_UNR[i,] = mvNormal(mod(betamat_unr, G[i,causal_ind]), cov_unr)
    return(P0_UNR)
    end
end


## Get Y1 and Y2 from phenotype matrix
function splitPhenotypeMatrix(P,nassoc1,nassoc2, npheno1,npheno2)
    n = size(P, 1)
    m = size(P, 2)
    if nassoc1>npheno1
        error("Error: nassoc1>npheno1")
    end
    if nassoc2>npheno2
        error("Error: nassoc2>npheno2")
    end
    if npheno1+npheno2 != m
        error("Error: npheno1+npheno2 != npheno")
    end
    if nassoc1 == 0
        Y2 = P[1:size(P,1), 1:npheno2] 
        Y1 = P[1:size(P,1), (npheno2+1):m]
    elseif nassoc2 == 0
        Y1 = P[1:size(P, 1), 1:npheno1]
        Y2 = P[1:size(P, 1),(npheno1+1):m]
    else 
        Y1 = P[:, vcat(1:nassoc1,(nassoc1+nassoc2+1):(nassoc1+nassoc2+(npheno1-nassoc1)))]
        Y2 = P[:, vcat((nassoc1+1):(nassoc1+nassoc2),(nassoc1+nassoc2+npheno1-nassoc1+1):m)]
    end
    return([Y1=Y1,Y2=Y2])
end


## Each row in Y1 is simulated as normal with the row of Y2 as mean (assuming beta=1), and cov=I
function simulateFullyMediatedPhenotype(Y2,npheno1,n_unrelated,traitcor)
    Y1 = zeros(npheno1*n_unrelated, npheno1)
    cov = createCovMatrix(npheno1, traitcor)
    for i in 1:n_unrelated
        if size(Y2, 2) <= npheno1
            mu = repeat(0,outer=npheno1)
            for j in 1:size(Y2, 2)
                mu[j] = Y2[i,j]
            end
        else 
            mu = Y2[i,1:npheno1]
        end
        Y1[i,] = MvNormal(mu, cov)
    end
    return(Y1)
end


function simulatePartiallyMediatedPhenotype(Y2,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno1, traitcor, causal_ind, nassoc1, variant, MAF_unr)
    betamat_unr = out2[:betamat]
    cov_unr = out2[:cov]
    ## Actual phenotype simulation:
    Y1 = zeros(npheno1*n_unrelated, npheno1)
    for i in 1:n_unrelated
        if size(Y2, 2) <= npheno1
            mu = repeat(0, outer=npheno1)
            for j in 1:size(Y2, 2)
                mu[j] = Y2[i, j]
            end
        else
            mu = Y2[i, 1:npheno1]
        end
        Y1[i,] = MvNormal(mod(betamat_unr, UNR_OBS[i,causal_ind]+mu), cov_unr)
    end
    return(Y1)
end


function simulatePhenotypesMediation(npheno1, npheno2,traitcor, nassoc1, nassoc2, causal_ind, MAF_unr, n_unrelated, variant,UNR_OBS, approach=1, numpcs=1, ignoreZ=false) 
    if approach == 1
            P0_UNR = simulatePhenotypes(npheno1+npheno2, traitcor, causal_ind, nassoc1+nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
            ## Now, we want to split this matrix into two to test for mediation
            out3 = splitPhenotypeMatrix(P0_UNR,nassoc1, nassoc2,npheno1,npheno2)
            Y1 = out3[:Y1]
            Y2 = out3[:Y2]
    elseif approach == 20 #No Mediatio 
        Y1 = simulatePhenotypes(npheno1, traitcor, causal_ind, nassoc1, variant, MAF_unr, n_unrelated, UNR_OBS)
        Y2 = simulatePhenotypes(npheno2, traitcor, causal_ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
    elseif approach == 2 && nassoc1 == 0 #Full Mediation 
        Y2 = simulatePhenotypes(npheno2, traitcor, causal_ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
        Y1 = simulateFullyMediatedPhenotype(Y2, npheno1, n_unrelated, traitcor)
    elseif approach == 2 && nassoc1 > 0 # Partial Mediation
        Y2 = simulatePhenotypes(npheno2, traitcor, causal_ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS)
        Y1 = simulatePartiallyMediatedPhenotype(Y2,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS)
    elseif approach == 3
        Z = [rnorm(n_unrelated)]
        if nassoc2 == 0 ##approach3: no association with G
            if nassoc1 != 0
                error("nassoc1 cannot be != 0 if nassoc2 is ==0")
            end
            P2 = simulateFullyMediatedPhenotype(Z, npheno2, n_unrelated, traitcor)
            P1 = simulateFullyMediatedPhenotype(Z, npheno1, n_unrelated, traitcor)
        elseif nassoc2 > 0
            P2 = simulatePartiallyMediatedPhenotype(Z,npheno2,n_unrelated, traitcor,causal_ind,nassoc2,variant,MAF_unr,UNR_OBS)
            if nassoc1 == 0
                P1 = simulateFullyMediatedPhenotype(Z, npheno1, n_unrelated, traitcor)
            elseif nassoc1 > 0
                P1 = simulatePartiallyMediatedPhenotype(Z,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS) 
            else
                error("nassoc1 should be >=0")
            end
        end
        if ignoreZ
            Y1 = P1
            Y2 = P2
        else
            ## PCA ---------------------------------------
            P = [P1 P2]
            #Distributions pkg
            P = scale(P)
            #MutlivariateStats pkg
            pc = pcasvd(P)
            prop = pc[sdev[1]]/sum(pc[sdev])
            print(["Proportion of variance explained by PC1: " prop])
            Z = [pc[x[1:size(x,1), 1:numpcs]]]
            ## ------------------------------------------- 
            Y2 = [convert(Int64, (size(P2, 1)*size(P2, 2))), size(P2, 2)]
            for i in 1:size(P2, 2)
                #GLM
                f = lm(P1[1, i]~Z)
##------------------------------------------------------------------------------
## Convert: 
##------------------------------------------------------------------------------
                Y2[:i] = residuals(f)
            end
            Y1 = [convert(Int64, (size(P1, 1)*size(P1, 2))), size(P1, 2)]
            for i in 1:size(P1, 2)
                f = lm(P1[1, i]~Z)
##------------------------------------------------------------------------------
## Convert: 
##------------------------------------------------------------------------------
                Y1[:i] = residuals(f)
            end
        end
    end
    return([Y1=Y1, Y2=Y2])
end


##------------------------------------------------------------------------------
## Simulating genotypes:
##------------------------------------------------------------------------------

UNR_OBS = [size(n_unrelated,1),repeat(NaN,outer=n_variants*n_unrelated)]
for i in 1:n_variants
    out = minorAlleleCountsBinomial(n_unrelated, maf)
    UNR_OBS[:i]= [out[G_mom]]
end
n_causal = floor(causal_var*n_variants)
causal_ind = sample(1:n_variants,n_causal,replace=false)

## determine the MAF of each variant in the sample
MAF_unr = mean(UNR_OBS, 2)/2

## we set "variant" to use the same functions as with cosi
if maf > 0.05
    variant = "common"
else
    variant = "rare"
end


##------------------------------------------------------------------------------
## Simulating phenotypes:
##------------------------------------------------------------------------------
out2 = simulatePhenotypesMediation(npheno1, npheno2,traitcor, nassoc1, nassoc2, causal_ind, MAF_unr, n_unrelated, variant, UNR_OBS,approach=sim_approach, ignoreZ=ignoreZ)
Y1 = out2[Y1]
Y2 = out2[Y2]
