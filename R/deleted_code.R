#' timed_epistasis_lik
#'
#' For a pair of variants (or two groups of variants), creates a likelihood function for a
#' model of pairwise epistasis with a "strong mutation, weak selection" assumption and an estimate of which mutation occurred first. 
#' 
#' The arguments to this function are automatically supplied by \code{ces_epistasis()} and \code{ces_gene_epistasis()}.
#' 
#' @param with_just_1 two-item list of baseline rates in v1/v2 for tumors with mutation in just the first variant(s)
#' @param with_just_2 two-item list of baseline rates in v1/v2 for tumors with mutation in just the second variant(s)
#' #' @param with_both two-item list of baseline rates in v1/v2 for tumors with mutation in both
#' @param with_neither two-item list of baseline rates in v1/v2 for tumors with mutation n neither
#' @param which_first three-item list for samples with mutation in both, specifying whether v1/v2/unspecified occurred first
#' 
#'
#' @export
#' @return A likelihood function


timed_epistasis_lik <- function(with_just_1, with_just_2, with_both, with_neither, which_first) 
{
  fn = function(par) {
    # sometimes the pars end up as NaNs or NAs, possibly because of inappropriate optimization techniques
    if(! all(is.finite(par))) {
      return(1e200)
    }
    
    # two points of discontinuity we need to account for
    if((par[3] == par[1] + par[2]) |
       (par[4] == par[1] + par[2])){return(1e200)}
    
    sum_log_lik <- 0
    
    if(! is.null(with_neither)) {
      # log(P{wt}) = -(A + B)
      A = par[1] * with_neither[[1]]
      B = par[2] * with_neither[[2]]
      ll = -1 * (A + B)
      sum_log_lik = sum_log_lik + sum(ll)  
    }
    
    
    if(! is.null(with_just_1)) {
      A = par[1] * with_just_1[[1]]
      B = par[2] * with_just_1[[2]]
      B_on_A = par[4] * with_just_1[[2]] 
      
      lik  = (A / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))
      sum_log_lik = sum_log_lik + sum(log(lik))   
    }
    
    
    if(! is.null(with_just_2)) {
      A = par[1] * with_just_2[[1]]
      B = par[2] * with_just_2[[2]]
      A_on_B = par[3] * with_just_2[[1]]
      
      lik = (B / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))
      sum_log_lik = sum_log_lik + sum(log(lik))
    }
    
    if(! is.null(with_both)) {
      
      if(length(which_first$v1 > 0)) {
        
        with_both_1b2 <- with_both[samples %in% which_first$v1]
        
        A = par[1] * with_both_1b2[[1]]
        B = par[2] * with_both_1b2[[2]]
        B_on_A = par[4] * with_both_1b2[[2]] 
        
        lik  = (A / (A + B)) * (1 - exp(-1 * (A + B)) - (((A + B) / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))))
        sum_log_lik = sum_log_lik + sum(log(lik))
      }
      
      if(length(which_first$v2 > 0)) {
        
        with_both_1b2 <- with_both[samples %in% which_first$v1]
        
        A = par[1] * with_both_2b1[[1]]
        B = par[2] * with_both_2b1[[2]]
        A_on_B = par[4] * with_both_2b1[[2]] 
        
        lik  = (B / (A + B)) * (1 - exp(-1 * (A + B)) - ( ((A + B) / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))))
        sum_log_lik = sum_log_lik + sum(log(lik))
      }
      
      if(length(which_first$unspecified > 0)) {
        
        with_both_unspecified <- with_both[samples %in% which_first$unspecified]
        
        A = par[1] * with_both_unspecified[[1]]
        B = par[2] * with_both_unspecified[[2]]
        A_on_B = par[3] * with_both_unspecified[[1]]
        B_on_A = par[4] * with_both_unspecified[[2]] 
        
        
        p_wt = exp(-1 * (A+B))
        p_A = (A / (A + B - B_on_A)) * (exp(-1 * B_on_A) - exp(-1 * (A + B)))
        p_B = (B / (A + B - A_on_B)) * (exp(-1 * A_on_B) - exp(-1 * (A + B)))   
        p_AB_unspecified = 1 - p_wt - p_A - p_B
        sum_log_lik = sum_log_lik + sum(log(p_AB_unspecified))
      }
    }
    
    
    # in case it tried all the max at once.
    if(!is.finite(sum_log_lik)){
      return(1e200)
    }
    return(-sum_log_lik)
  }
  
  # Set default values for all parameters, which ces_variant will use to set starting values of optimization
  formals(fn)[["par"]] = 1000:1003
  
  # Optimization tool, bbmle::mle, requires that vector of parameters to optimize have named elements
  bbmle::parnames(fn) = c("ces_v1", "ces_v2", "ces_v1_after_v2", "ces_v2_after_v1")
  return(fn)
}