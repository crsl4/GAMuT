##------------------------------------------------------------------------------
## Variables:  
##------------------------------------------------------------------------------
maf = 0.25 
traitcor = "med"
nassoc1 = 1 
nassoc2 = 1 
npheno1 = 1 
npheno2 = 1 
n_unrelated = 100 
n_variants = 500
causal_var = 0.01
test_approach = 1 
sim_approach = 3 
ignoreZ = true

npheno = npheno1 + npheno2
nassoc = nassoc1
variant = "common"

#UNR_OBS = fill(NaN, 100, 500*100)
#for i in 1:n_variants
    #out = minorAlleleCountsBinomial(n_unrelated, maf)
    #UNR_OBS[:i] = [out[:eur_mom]]
#end

#MAF_unr = sum(UNR_OBS, 2)/2
MAF_unr = rand(500,1)

n_causal = floor(causal_var*n_variants)
n_causal = 5
using StatsBase
causal_ind = sample(1:n_variants,n_causal,replace=false)
causal_ind = [407 448 142 340 261]
