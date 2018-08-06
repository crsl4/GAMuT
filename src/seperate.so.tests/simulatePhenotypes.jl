function simulatePhenotypes(npheno, traitcor, causal_ind, nassoc, variant, MAF_unr, n_unrelated, G, effectSize)
    ## Creating parameters to simulate phenotype
    out2 = parameters4phenotypeSimulation(npheno, traitcor, causal_ind, nassoc, variant, MAF_unr, effectSize)
    betamat_unr = out2[1]
    cov_unr = out2[2]
    ## Actual phenotype simulation:
    P0_UNR = zeros(n_unrelated,npheno)
    for i in 1:n_unrelated
        mu = betamat_unr * G[i,causal_ind]
        mu = vcat(mu...) ## splat converts array to iterable
        P0_UNR[i,:] = rand(MvNormal(mu, cov_unr),1)
    end
    return(P0_UNR)
end


