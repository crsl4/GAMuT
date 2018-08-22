## This function is not used in our simulation. We simply use simulatePhenotypes.
## In the PCA part of this function, the scale and pcasvd functions are deprecated. pcasvd is part of the MutlivariateStats package which is deprecated.
## The lm function is also deprecated due to the deprecation of the GLM package. 
## This function has not been fully tested due to the deprecation of the above functions, and it does not fully run. 

using GLM
using MultivariateStats 
function simulatePhenotypesMediation(npheno1, npheno2,traitcor, nassoc1, nassoc2, causal_ind, MAF_unr, n_unrelated, variant,UNR_OBS, effectSize, approach, numpcs, ignoreZ) ##approach and numpcs = 1, ignoreZ = false
    if approach == 1
            P0_UNR = simulatePhenotypes(npheno1+npheno2, traitcor, causal_ind, nassoc1+nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS, effectSize)
            ## Now, we want to split this matrix into two to test for mediation
            out3 = splitPhenotypeMatrix(P0_UNR,nassoc1, nassoc2,npheno1,npheno2)
            Y1 = out3[:Y1]
            Y2 = out3[:Y2]
    elseif approach == 20 #No Mediatio
        Y1 = simulatePhenotypes(npheno1, traitcor, causal_ind, nassoc1, variant, MAF_unr, n_unrelated, UNR_OBS, effectSize)
        Y2 = simulatePhenotypes(npheno2, traitcor, causal_ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS, effectSize)
    elseif approach == 2 && nassoc1 == 0 #Full Mediation
        Y2 = simulatePhenotypes(npheno2, traitcor, causal_ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS, effectSize)
        Y1 = simulateFullyMediatedPhenotype(Y2, npheno1, n_unrelated, traitcor)
    elseif approach == 2 && nassoc1 > 0 # Partial Mediation
        Y2 = simulatePhenotypes(npheno2, traitcor, causal_ind, nassoc2, variant, MAF_unr, n_unrelated, UNR_OBS, effectSize)
        Y1 = simulatePartiallyMediatedPhenotype(Y2,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS,effectSize)
    elseif approach == 3
        Z = [rnorm(n_unrelated)]
        if nassoc2 == 0 ##approach3: no association with G
            if nassoc1 != 0
                error("nassoc1 cannot be != 0 if nassoc2 is ==0")
            end
            P2 = simulateFullyMediatedPhenotype(Z, npheno2, n_unrelated, traitcor)
            P1 = simulateFullyMediatedPhenotype(Z, npheno1, n_unrelated, traitcor)
        elseif nassoc2 > 0
            P2 = simulatePartiallyMediatedPhenotype(Z,npheno2,n_unrelated, traitcor,causal_ind,nassoc2,variant,MAF_unr,UNR_OBS,effectSize)
            if nassoc1 == 0
                P1 = simulateFullyMediatedPhenotype(Z, npheno1, n_unrelated, traitcor)
            elseif nassoc1 > 0
                P1 = simulatePartiallyMediatedPhenotype(Z,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS,effectSize)
            else
                error("nassoc1 should be >=0")
            end
        end
        if ignoreZ
            Y1 = P1
            Y2 = P2
        else
            ## PCA ---------------------------------------
            P = hcat(P1, P2)
            #Distributions pkg
            P = scale(P)
            #MutlivariateStats pkg
            pc = pcasvd(P)
            prop = pc[:sdev[1]]/sum(pc[:sdev])
            print(["Proportion of variance explained by PC1: " prop])
            Z = [pc[:x[1:size(:x,1), 1:numpcs]]]
            ## -------------------------------------------
            Y2 = [convert(Int64, (size(P2, 1)*size(P2, 2))), size(P2, 2)]
            for i in 1:size(P2, 2)
                #GLM
                f = lm(P1[1, i]~Z)
                Y2[:i] = StatsBase.residuals(f)
            end
            Y1 = [convert(Int64, (size(P1, 1)*size(P1, 2))), size(P1, 2)]
            for i in 1:size(P1, 2)
                f = lm(P1[1, i]~Z)
                Y1[:i] = StatsBase.residuals(f)
            end
        end
    end
    return([Y1, Y2])
end
