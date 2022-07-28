#' sswm_lik
#'
#' Generates log-likelihood function of site-level selection with "strong selection, weak
#' mutation" assumption. All arguments to this likelihood function factory are
#' automatically supplied by \code{ces_variant()}.
#'
#' @param rates_tumors_with vector of site-specific mutation rates for all tumors with variant
#' @param rates_tumors_without vector of site-specific mutation rates for all eligible tumors without variant
#' @export
sswm_lik = function(rates_tumors_with, rates_tumors_without) {
  fn = function(gamma) {
    gamma = unname(gamma) # math faster on unnamed vectors
    sum_log_lik = 0
    if (length(rates_tumors_without) > 0) {
      sum_log_lik = -1 * sum(gamma * rates_tumors_without)
    }
    if (length(rates_tumors_with) > 0) {
      sum_log_lik = sum_log_lik + sum(log(1 - exp(-1 * gamma * rates_tumors_with)))
    }

    # convert to negative loglikelihood and return
    return(-1 * sum_log_lik)
  }
  
  # Set default values for gamma (SI), which ces_variant will use to set starting value of optimization
  formals(fn)[["gamma"]] = 1000
  bbmle::parnames(fn) = "selection_intensity"
  return(fn)
}
  
  
#' sswm_sequential_lik
#' 
#' As in sswm_lik, selection intensities are calculated at variant sites under a "strong
#' selection, weak mutation" assumption. In this version, each sample is assigned to one
#' of an ordered set of disease progression states, and selection is assumed to vary
#' across states. For example, in a two-state local/metastatic model, each variant has two
#' independent selection intensities. Metastatic samples could have acquired the variant
#' while in their current state or at some earlier time, while the local state selection
#' intensity applied.
#' 
#' All arguments to this likelihood function factory are automatically supplied by
#' \code{ces_variant()}.
#' 
#' @param rates_tumors_with named vector of site-specific mutation rates for all tumors
#'   with variant
#' @param rates_tumors_without named vector of site-specific mutation rates for all
#'   eligible tumors without variant
#' @param sample_index data.table with columns Unique_Patient_Identifier, group_name, group_index
#' @keywords internal
sswm_sequential_lik <- function(rates_tumors_with, rates_tumors_without, sample_index) {
  stages_tumors_with = sample_index[names(rates_tumors_with), group_index, on = "Unique_Patient_Identifier"]
  stages_tumors_without = sample_index[names(rates_tumors_without), group_index, on = "Unique_Patient_Identifier"]

  # e.g., 1 -> 1; 3 -> 1,2,3 (for vector subsetting)
  sequential_stages_tumors_with = list() # unused if no tumors
  if(length(stages_tumors_with) > 0) {
    sequential_stages_tumors_with = lapply(stages_tumors_with, function(x) 1:x) 
  }
  num_pars = sample_index[, uniqueN(group_index)]
  
  fn = function(gamma) {
    gamma = unname(gamma)
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
    
    if(length(rates_tumors_with) > 0) {
      gamma_sums = mapply(calc_gamma_sums_mut, rates_tumors_with, sequential_stages_tumors_with)
      sum_log_lik = sum_log_lik + sum(gamma_sums)
    }
    
    # in case it tried all the max at once
    if(!is.finite(sum_log_lik)){
      return(-1e200)
    }
    return(-1 * sum_log_lik)
  }

  # Set default values for all parameters, which ces_variant will use to set starting values of optimization
  formals(fn)[["gamma"]] = rep.int(1000, num_pars)
  
  # Optimization tool, bbmle::mle, requires that vector of parameters to optimize have named elements
  group_names = unique(sample_index[, .(group_index, group_name)], by = 'group_index')[order(group_index), group_name]
  bbmle::parnames(fn) = paste0("si_", group_names)
  return(fn)
}


