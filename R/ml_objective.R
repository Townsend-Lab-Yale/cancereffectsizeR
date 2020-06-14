#' ml_objective
#'
#' Log-likelihood function under a model of site-specific selection intensities at variant sites
#'
#' @param gamma A selection intensity at which to calculate the likelihood
#' @param tumor_stages an environment with keys = tumor names, values = stage of tumor
#' @param tumors_without_gene_mutated vector of (eligible) tumors without any mutation in the gene of the variant
#' @param tumors_with_variant vector of (eligible) tumors with the variant
#' @param baseline_rates a named vector of baseline mutation rates for each eligible tumor
#' @return A log likelihood value
#' @export
#' @keywords internal

ml_objective <- function(gamma, tumor_stages, tumors_with_variant, tumors_without_gene_mutated, baseline_rates, modifier = 0) {
  
  fn = function(gamma) {
    sums = cumsum(gamma)
    gamma_sums = sums[tumor_stages[tumors_without_gene_mutated]]
    sum_log_lik = -1 * sum(gamma_sums * baseline_rates[tumors_without_gene_mutated])

    calc_gamma_sums_mut = function(tumor) {
      # stage-specific likelihoods of mutation
      lik_no_mutation = exp(-1 * gamma * baseline_rates[tumor])
      lik_mutation = 1 - lik_no_mutation
      
      this_sum = lik_mutation[1] # likelihood of mutation having occurred by stage 1 
      
      # for tumors that are at stage > 1, add in likelihoods of mutation having occurred at later stages
      current_stage = 2
      while (current_stage <= tumor_stages[[tumor]]) {
        # e.g., at current_stage = 3, take product of no mutation in stage 1, no mutation in stage 2, yes mutation in stage 3
        this_sum = this_sum + lik_mutation[current_stage] * prod(lik_no_mutation[(current_stage-1):1])
        current_stage <- current_stage + 1
      }
      return(log(this_sum))
    }
    
    gamma_sums = sapply(tumors_with_variant, calc_gamma_sums_mut)
    sum_log_lik = sum_log_lik + sum(gamma_sums)
    
    # in case it tried all the max at once.
    if(!is.finite(sum_log_lik)){
      return(-1e200)
    }
    return(-1 * sum_log_lik - modifier)
  }
  
  return(fn)
    
}
