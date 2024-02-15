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
  formals(fn)[["gamma"]] = 1
  bbmle::parnames(fn) = "selection_intensity"
  return(fn)
}
  


