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
ml_objective_epistasis_full_gene <- function(par,
                                             MAF_input1 = NULL,
                                             MAF_input2,
                                             all_tumors,
                                             gene1,
                                             gene2,
                                             specific_mut_rates1,
                                             specific_mut_rates2,
                                             variant_freq_1,
                                             variant_freq_2) {

  # if optimization algorithm is not calling this function correctly, arguments end up null
  if (is.null(MAF_input1)) {
    stop("ml_objective_epistasis_full_gene called without required arguments")
  }
  
  # sometimes the pars end up as NaNs or NAs, possibly because of inappropriate optimization techniques
  if(! all(is.finite(par))) {
    return(1e200)
  }
  # only the tumors containing a recurrent variant factor into the selection analysis
  tumors_with_variant1_mutated <- unique(MAF_input1[which(MAF_input1$unique_variant_ID_AA %in% names(which(variant_freq_1>1))),"Unique_Patient_Identifier"])
  tumors_with_variant2_mutated <- unique(MAF_input2[which(MAF_input2$unique_variant_ID_AA %in% names(which(variant_freq_2>1))),"Unique_Patient_Identifier"])

  tumors_with_both_mutated <- base::intersect(tumors_with_variant1_mutated,tumors_with_variant2_mutated)

  tumors_with_ONLY_variant1_mutated <- tumors_with_variant1_mutated[which(!tumors_with_variant1_mutated %in% tumors_with_variant2_mutated)]

  tumors_with_ONLY_variant2_mutated <- tumors_with_variant2_mutated[which(!tumors_with_variant2_mutated %in% tumors_with_variant1_mutated)]

  tumors_with_neither_mutated <- setdiff(all_tumors, c(tumors_with_both_mutated,tumors_with_ONLY_variant1_mutated,tumors_with_ONLY_variant2_mutated))



  # two points of discontinuity we need to account for
  if((par[3] == par[1] + par[2]) |
     (par[4] == par[1] + par[2])){return(1e200)}

  sum_log_lik <- 0

  # not in either

  for (tumor in tumors_with_neither_mutated) {
    sum_log_lik <- sum_log_lik +
      (
        (-1*(par[1] * sum(specific_mut_rates1[tumor, ]))) +
          (-1*(par[2] * sum(specific_mut_rates2[tumor, ])))
      )

  }


  # tumors only have variant 1


  for( tumor in tumors_with_ONLY_variant1_mutated){
    sum_log_lik <- sum_log_lik +

      # if(
      log(
        # log(
        -1*(
          (par[1] * sum(specific_mut_rates1[tumor, ])) /
            (
              (par[1] * sum(specific_mut_rates1[tumor, ])) +
                (par[2] * sum(specific_mut_rates2[tumor, ])) -
                (par[4] * sum(specific_mut_rates2[tumor, ]))
            )
          # )
        ) *
          # log(
          ((exp(  (-1*(par[1] * sum(specific_mut_rates1[tumor, ]))) +
                    (-1*(par[2] * sum(specific_mut_rates2[tumor, ]))))) -
             (exp(-1*(par[4] * sum(specific_mut_rates2[tumor, ])))))
        # )
      )

  }

  # tumors only have variant 2

  for( tumor in tumors_with_ONLY_variant2_mutated){
    sum_log_lik <- sum_log_lik +
      log(
        # log(
        -1* (
          (par[2] * sum(specific_mut_rates2[tumor, ])) /
            (
              (par[2] * sum(specific_mut_rates2[tumor, ])) +
                (par[1] * sum(specific_mut_rates1[tumor, ])) -
                (par[3] * sum(specific_mut_rates1[tumor, ]))
            )
          # )
        ) *
          # log(
          ((exp(  (-1*(par[1] * sum(specific_mut_rates1[tumor, ]))) +
                    (-1*(par[2] * sum(specific_mut_rates2[tumor, ]))))) -
             (exp(-1*(par[3] * sum(specific_mut_rates1[tumor, ])))))
        # )
      )

  }


  # tumors have both variants
  for( tumor in tumors_with_both_mutated) {

    sum_log_lik <- sum_log_lik + log(
      1-
        (
          ( # P(wt)
            exp((-1*(par[1] * sum(specific_mut_rates1[tumor, ]))) +
                  (-1*(par[2] * sum(specific_mut_rates2[tumor, ]))))
          ) +
            ( #P(1)
              # log(
              -1*(
                (par[1] * sum(specific_mut_rates1[tumor, ])) /
                  (
                    (par[1] * sum(specific_mut_rates1[tumor, ])) +
                      (par[2] * sum(specific_mut_rates2[tumor, ])) -
                      (par[4] * sum(specific_mut_rates2[tumor, ]))
                  )
                # )
              ) *
                # log(
                ((exp(  (-1*(par[1] * sum(specific_mut_rates1[tumor, ]))) +
                          (-1*(par[2] * sum(specific_mut_rates2[tumor, ]))))) -
                   (exp(-1*(par[4] * sum(specific_mut_rates2[tumor, ])))))
              # )
            ) +
            ( # P(2)
              -1* (
                (par[2] * sum(specific_mut_rates2[tumor, ])) /
                  (
                    (par[2] * sum(specific_mut_rates2[tumor, ])) +
                      (par[1] * sum(specific_mut_rates1[tumor, ])) -
                      (par[3] * sum(specific_mut_rates1[tumor, ]))
                  )
                # )
              ) *
                # log(
                ((exp(  (-1*(par[1] * sum(specific_mut_rates1[tumor, ]))) +
                          (-1*(par[2] * sum(specific_mut_rates2[tumor, ]))))) -
                   (exp(-1*(par[3] * sum(specific_mut_rates1[tumor, ])))))
              # )
            )


        )


    )



  }


  # in case it tried all the max at once.
  if(!is.finite(sum_log_lik)){
    return(1e200)
  }
  return(-sum_log_lik)
}
