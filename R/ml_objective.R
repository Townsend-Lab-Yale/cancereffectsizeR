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
  
  stages_tumors_without = unname(tumor_stages[tumors_without_gene_mutated])
  rates_tumors_without = unname(baseline_rates[tumors_without_gene_mutated])
  
  stages_tumors_with = unname(tumor_stages[tumors_with_variant])
  sequential_stages_tumors_with = lapply(stages_tumors_with, function(x) 1:x) # e.g., 1 -> 1; 3 -> 1,2,3 (for vector subsetting)
  rates_tumors_with = unname(baseline_rates[tumors_with_variant])
  
  fn = function(gamma) {
    sums = cumsum(gamma)
    gamma_sums = sums[stages_tumors_without]
    sum_log_lik = -1 * sum(gamma_sums * rates_tumors_without)

    calc_gamma_sums_mut = function(rate, stage_seq) {
      # stage-specific likelihoods of mutation
      lik_no_mutation = exp(-1 * gamma * rate)
      lik_mutation = 1 - lik_no_mutation
      
      cum_lik_no_mut = c(1, cumprod(lik_no_mutation))
      return(log(sum(cum_lik_no_mut[stage_seq] * lik_mutation[stage_seq])))
    }
    
    if(length(tumors_with_variant) > 0) {
      gamma_sums = mapply(calc_gamma_sums_mut, rates_tumors_with, sequential_stages_tumors_with)
      sum_log_lik = sum_log_lik + sum(gamma_sums)
    }
    
    # in case it tried all the max at once.
    if(!is.finite(sum_log_lik)){
      return(-1e200)
    }
    return(-1 * sum_log_lik - modifier)
  }
  
  return(fn)
    
}
