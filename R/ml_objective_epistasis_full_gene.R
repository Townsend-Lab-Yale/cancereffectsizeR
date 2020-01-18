#' ml_objective for epistasis calculation at the gene level
#'
#' Objective function that we will be optimizing (minimizing the negative of the log likelihood, so maximizing the log likelihood) in order to find the site specific selection intensity that maximizes the likelihood of each tumor having a mutation or not, where the mutation rates are site and tumor specific.
#'
#'
#' @param par A selection intensity at which to calculate the likelihood
#' @param MAF_input A data frame that includes columns "Unique_patient_identifier", "Gene_name", and "unique_variant_ID_AA"
#' @param all_tumors A list of all the tumors we are calculating the likelihood across
#' @param gene The gene we want to look at
#' @param variant The variant we want to look at
#' @param specific_mut_rates A matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants (produced by mutation_rate_calc)
#'
#' @return A log likelihood value
#' @export
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
  
  tmp = with_neither[,(-1*par[1]*dt1_sum) + (-1 * par[2] * dt2_sum) ]
  sum_log_lik = sum_log_lik + sum(tmp)
  
  tmp = with_just_1[,log((-1 * (par[1] * dt1_sum / (par[1]*dt1_sum + par[2]*dt2_sum - par[4]*dt2_sum))) * (exp((-1 * par[1]*dt1_sum) + (-1 * par[2]* dt2_sum)) - exp(-1*par[4] * dt2_sum)))]
  sum_log_lik = sum_log_lik + sum(tmp)
  
  tmp = with_just_2[,log((-1 * (par[2] * dt2_sum / (par[2]*dt2_sum + par[1]*dt1_sum - par[3]*dt1_sum))) * (exp((-1 * par[1]*dt1_sum) + (-1 * par[2]* dt2_sum)) - exp(-1*par[3] * dt1_sum)))]
  sum_log_lik = sum_log_lik + sum(tmp)
  tmp = with_both[,
                  log(
                    1-
                      (
                        ( # P(wt)
                          exp((-1*(par[1] * dt1_sum)) +
                                (-1*(par[2] * dt2_sum)))
                        ) +
                          ( #P(1)
                            # log(
                            -1*(
                              (par[1] * dt1_sum) /
                                (
                                  (par[1] * dt1_sum) +
                                    (par[2] * dt2_sum) -
                                    (par[4] * dt2_sum)
                                )
                              # )
                            ) *
                              # log(
                              ((exp(  (-1*(par[1] * dt1_sum)) +
                                        (-1*(par[2] * dt2_sum)))) -
                                 (exp(-1*(par[4] * dt2_sum))))
                            # )
                          ) +
                          ( # P(2)
                            -1* (
                              (par[2] * dt2_sum) /
                                (
                                  (par[2] * dt2_sum) +
                                    (par[1] * dt1_sum) -
                                    (par[3] * dt1_sum)
                                )
                              # )
                            ) *
                              # log(
                              ((exp(  (-1*(par[1] * dt1_sum)) +
                                        (-1*(par[2] * dt2_sum)))) -
                                 (exp(-1*(par[3] * dt1_sum))))
                            # )
                          )
                        
                        
                      )
                  )
                  ]
  sum_log_lik = sum_log_lik + sum(tmp)


  # in case it tried all the max at once.
  if(!is.finite(sum_log_lik)){
    return(1e200)
  }
  return(-sum_log_lik)
}
