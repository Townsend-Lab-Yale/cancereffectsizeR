#' Calculate selection intensity and epistasis
#'
#' Actual function to find the site specific selection intensity that maximizes the likelihood of each tumor being mutated or not. Uses site and tumor specific mutation rates. Uses Brent 1 dimensional optimization technique.optimize_gamma
#'
#' @param MAF_input A data frame that includes columns "Unique_patient_identifier", "Gene_name", and "unique_variant_ID_AA"
#' @param all_tumors A list of all the tumors we are calculating the likelihood across
#' @param gene The gene we want to look at
#' @param variant The variant we want to look at
#' @param specific_mut_rates A matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants (produced by mutation_rate_calc)
#'
#' @return The optimal selection intensity for the gene and variant that maximizes the likelihood of the observations
#' @export
#'
#' @examples
optimize_gamma_epistasis <- function(MAF_input1,
                           MAF_input2,
                           all_tumors,
                           gene1,
                           gene2,
                           variant1,
                           variant2,
                           specific_mut_rates1,
                           specific_mut_rates2) {

  # par_init <- rep(1000,length=4)
  par_init <- 1000:1003

  return(optim(par=par_init,
               fn=cancereffectsizeR::ml_objective_epistasis,
               MAF_input1=MAF_input1,
               MAF_input2=MAF_input2,
               all_tumors=all_tumors,
               gene1=gene1,
               gene2=gene2,
               variant1=variant1,
               variant2=variant2,
               specific_mut_rates1=specific_mut_rates1,
               specific_mut_rates2=specific_mut_rates2,
               method="L-BFGS-B",
               lower=1e-6,
               upper=1e7,
               control=list(fnscale=-1e-12))$par)
}
