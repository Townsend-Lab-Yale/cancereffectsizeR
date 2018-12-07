
# Actual function to find the site specific selection intensity that maximizes the likelihood of each tumor
# being mutated or not. Uses site and tumor specific mutation rates. Uses Brent 1 dimensional optimization technique.

# Inputs:
# all_tumors: a list of all the tumors we are calculating the likelihood across
# gene, variant: the gene and variant we want to look at
# specific_mut_rates: a matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants

# Outputs:
# The optimal selection intensity for the gene and variant that maximizes the likelihood of the observations

optimize_gamma <- function(MAF_input, all_tumors, gene, variant, specific_mut_rates) {
  return(optim(par=1000, fn=cancereffectsizeR::ml_objective, MAF_input=MAF_input, all_tumors=all_tumors, gene=gene, variant=variant, specific_mut_rates=specific_mut_rates,
               method="Brent", lower=1, upper=1000000000, control=list(fnscale=-1))$par)
}
