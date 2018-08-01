function simulatePartiallyMediatedPhenotype(Y2,npheno1,n_unrelated, traitcor,causal_ind,nassoc1,variant,MAF_unr,UNR_OBS, effectSize)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno1, traitcor, causal_ind, nassoc1, variant, MAF_unr, effectSize)
    betamat_unr = out2[1]
    cov_unr = out2[2]
    ## Actual phenotype simulation:
    Y1 = zeros(n_unrelated, npheno1)
    for i in 1:n_unrelated
        if size(Y2, 2) <= npheno1
            mu = zeros(npheno1)
            for j in 1:size(Y2,2)
                mu = Y2[i, j]
            end
        else
            #Works
            mu = Y2[i, 1:npheno1]
        end
        mu0 = betamat_unr * UNR_OBS[i,causal_ind]'
        mu0 = vcat(mu0...)
        Y1[i,:] = rand(MvNormal( mu0 + mu, cov_unr),1)
    end
    return(Y1)
end
