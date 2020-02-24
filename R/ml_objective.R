#' ml_objective
#'
#' Objective function that we will be optimizing in order to find the site specific selection intensity that maximizes the likelihood of each tumor having a mutation or not, where the mutation rates are site and tumor specific.
#'
#' @param gamma A selection intensity at which to calculate the likelihood
#' @param tumor_stages an environment with keys = tumor names, values = stage of tumor
#' @param variant the unique mutation whose selection intensity is being determined
#' @param tumors_without_gene_mutated list of tumors without any mutation in the gene of the variant
#' @param tumors_with_pos_mutated list of tumors with the specific variant in question
#' @param specific_mut_rates A matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants (produced by mutation_rate_calc)
#' @return A log likelihood value
#'
#' @examples
ml_objective <- function(gamma, tumor_stages, tumors_without_gene_mutated, tumors_with_pos_mutated, 
  variant, specific_mut_rates, modifier=0) {

  calc_gamma_sums_no_mut = function(tumor) {
    return(sum(-1*gamma[1:tumor_stages[tumor]]))
  }

  gamma_sums = vapply(tumors_without_gene_mutated, calc_gamma_sums_no_mut, FUN.VALUE = 1.0)
  sum_log_lik = sum(gamma_sums * specific_mut_rates[tumors_without_gene_mutated, variant])


  calc_gamma_sums_mut = function(tumor) {
    # stage-specific likelihoods of mutation
    lik_no_mutation = exp(-1 * gamma * specific_mut_rates[tumor, variant])
    lik_mutation = 1 - lik_no_mutation

    this_sum = lik_mutation[1] # likelihood of mutation having occurred by stage 1 

    # for tumors that are at stage > 1, add in likelihoods of mutation having occurred at later stages
    current_stage = 2
    while (current_stage <= tumor_stages[[tumor]]) {
      # e.g., at current_stage = 3, take product of no mutation in stage 1, no mutation in stage 2, yes mutation in stage 3
      this_sum = this_sum + lik_mutation[current_stage] * prod(lik_no_mutation[current_stage-1:1])
      current_stage <- current_stage + 1
    }
    return(log(this_sum))
  }

  gamma_sums = vapply(tumors_with_pos_mutated, calc_gamma_sums_mut, FUN.VALUE = 1.0)
  sum_log_lik = sum_log_lik + sum(gamma_sums)

  # in case it tried all the max at once.
  if(!is.finite(sum_log_lik)){
    return(-1e200)
  }
  return(sum_log_lik - modifier)
}
