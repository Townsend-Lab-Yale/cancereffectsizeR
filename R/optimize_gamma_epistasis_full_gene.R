#' Calculate selection intensity and epistasis at the gene level
#'
#' Actual function to find the site specific selection intensity that maximizes the likelihood of each tumor being mutated or not. Uses site and tumor specific mutation rates. Uses Brent 1 dimensional optimization technique.optimize_gamma
#'
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
optimize_gamma_epistasis_full_gene <- function(MAF_input1,
                                               MAF_input2,
                                               all_tumors,
                                               gene1,
                                               gene2,
                                               specific_mut_rates1,
                                               specific_mut_rates2,
                                               variant_freq_1,
                                               variant_freq_2,
                                               full_gene_epistasis_lower_optim,
                                               full_gene_epistasis_upper_optim,
                                               full_gene_epistasis_fnscale) {

  # par_init <- rep(1000,length=4)
  par_init <- 1000:1003

  # Some methods crash for unclear reasons: Rnmin, nmkb, newuoa
  # Others should be skipped automatically because they aren't appropriate, but it's
  # necessary to skip them manually: "subplex", "snewtonm", "snewton", "CG", "BFGS","Nelder-Mead", "nlm", "lbfgs"
  # Leaving out Rtnmin method because it crashes (still unclear why)

  # This list consists of all other methods offered by optimx; a couple are much slower than others and should be dropped unless they're good
  methods = c("L-BFGS-B", "nlminb", "lbfgsb3", "Rcgmin", "Rvmmin",
              "spg", "bobyqa", "hjkb", "hjn")

  # To test/debug add this option to optimx options: "control = list(trace = 1) "
  par = 1000:1003
  
  # suppress warnings because user doesn't need to be informed about methods that fail to converge in particular instances
  opm_output <- suppressWarnings(optimx::opm(par,
                            fn=cancereffectsizeR::ml_objective_epistasis_full_gene,
                            method=methods,
                            lower=1e-3,
                            upper=1e20, gr = "grfwd",
                            MAF_input1=MAF_input1,
                            MAF_input2=MAF_input2,
                            all_tumors=all_tumors,
                            gene1=gene1,
                            gene2=gene2,
                            specific_mut_rates1=specific_mut_rates1,
                            specific_mut_rates2=specific_mut_rates2,
                            variant_freq_1=variant_freq_1,
                            variant_freq_2=variant_freq_2))

  #TODO: explore the best possible optimization algorithm, fnscale, etc.
  opm_output$value <- -opm_output$value # we did a minimization of the negative, so need to take the negative again to find the maximum
  opm_output <- opm_output[order(opm_output$value,decreasing = T),]

  return(opm_output[1,1:4])
}
