#' ml_objective for epistasis calculation at the gene level
#'
#' Objective function that we will be optimizing (minimizing the negative of the log likelihood, so maximizing the log likelihood) in order to find the site specific selection intensity that maximizes the likelihood of each tumor having a mutation or not, where the mutation rates are site and tumor specific.
#'
#'
#' @param par A selection intensity at which to calculate the likelihood
#' @param with_just_1
#' @param with_just_2
#' @param with_both
#' @param with_neither
#'
#' @return A log likelihood value
#'
#' @examples
ml_objective_epistasis_full_gene <- function(par, with_just_1, with_just_2, with_both, with_neither) {

  
  # sometimes the pars end up as NaNs or NAs, possibly because of inappropriate optimization techniques
  if(! all(is.finite(par))) {
    return(1e200)
  }
 


  # two points of discontinuity we need to account for
  if((par[3] == par[1] + par[2]) |
     (par[4] == par[1] + par[2])){return(1e200)}

  sum_log_lik <- 0
  
  if(! is.null(with_just_1)) {
    tmp = (-1*par[1]*with_neither[[1]] + (-1 * par[2] * with_neither[[2]]))
    sum_log_lik = sum_log_lik + sum(tmp)  
  }
  
  
  if(! is.null(with_just_1)) {
    tmp = log((-1 * (par[1] * with_just_1[[1]] / (par[1]*with_just_1[[1]] + par[2]*with_just_1[[2]] - par[4]*with_just_1[[2]]))) * (exp((-1 * par[1]*with_just_1[[1]]) + (-1 * par[2]* with_just_1[[2]])) - exp(-1*par[4] * with_just_1[[2]])))
    sum_log_lik = sum_log_lik + sum(tmp)   
  }
  
  
  if(! is.null(with_just_2)) {
    tmp = log((-1 * (par[2] * with_just_2[[2]] / (par[2]*with_just_2[[2]] + par[1]*with_just_2[[1]] - par[3]*with_just_2[[1]]))) * (exp((-1 * par[1]*with_just_2[[1]]) + (-1 * par[2]* with_just_2[[2]])) - exp(-1*par[3] * with_just_2[[1]])))
    sum_log_lik = sum_log_lik + sum(tmp)   
  }
  
  
  if(! is.null(with_both)) {
    tmp =  log(
      1-
        (
          ( # P(wt)
            exp((-1*(par[1] * with_both[[1]])) +
                  (-1*(par[2] * with_both[[2]])))
          ) +
            ( #P(1)
              # log(
              -1*(
                (par[1] * with_both[[1]]) /
                  (
                    (par[1] * with_both[[1]]) +
                      (par[2] * with_both[[2]]) -
                      (par[4] * with_both[[2]])
                  )
                # )
              ) *
                # log(
                ((exp(  (-1*(par[1] * with_both[[1]])) +
                          (-1*(par[2] * with_both[[2]])))) -
                   (exp(-1*(par[4] * with_both[[2]]))))
              # )
            ) +
            ( # P(2)
              -1* (
                (par[2] * with_both[[2]]) /
                  (
                    (par[2] * with_both[[2]]) +
                      (par[1] * with_both[[1]]) -
                      (par[3] * with_both[[1]])
                  )
                # )
              ) *
                # log(
                ((exp(  (-1*(par[1] * with_both[[1]])) +
                          (-1*(par[2] * with_both[[2]])))) -
                   (exp(-1*(par[3] * with_both[[1]]))))
              # )
            )
          
          
        )
    )
    
    sum_log_lik = sum_log_lik + sum(tmp)
  }
  


  # in case it tried all the max at once.
  if(!is.finite(sum_log_lik)){
    return(1e200)
  }
  return(-sum_log_lik)
}
