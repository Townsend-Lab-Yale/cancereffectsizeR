#' Calculate selection intensity
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
optimize_gamma <- function(MAF_input, all_tumors, gene, variant, specific_mut_rates) {
  par_init <- rep(1000,length=length(unique(all_tumors[,"subset"])))
    return(optim(par=par_init, fn=cancereffectsizeR::ml_objective, MAF_input=MAF_input, all_tumors=all_tumors, gene=gene, variant=variant, specific_mut_rates=specific_mut_rates,
                 method="L-BFGS-B", lower=1e-3, upper=1000000000, control=list(fnscale=-1e-12)))

}
