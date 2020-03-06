#' ml_objective for epistasis calculation
#'
#' Objective function that we will be optimizing in order to find the site specific selection intensity that maximizes the likelihood of each tumor having a mutation or not, where the mutation rates are site and tumor specific.
#'
#' @param gamma A selection intensity at which to calculate the likelihood
#' @param MAF_input1 A data frame that includes columns "Unique_patient_identifier", "Gene_name", and "unique_variant_ID_AA"
#' @param MAF_input2 A data frame that includes columns "Unique_patient_identifier", "Gene_name", and "unique_variant_ID_AA"
#' @param all_tumors A list of all the tumors we are calculating the likelihood across
#' @param variant1 one variant in epistatic analysis
#' @param variant2 the other
#' @param specific_mut_rates1 site and tumor specific mutation rates (matrix: rows = tumors, columns = variants; produced by mutation_rate_calc)
#' @param specific_mut_rates2 site and tumor specific mutation rates (matrix: rows = tumors, columns = variants; produced by mutation_rate_calc)
#'
#' @return A log likelihood value
#'
ml_objective_epistasis <- function(gamma, MAF_input1,MAF_input2, all_tumors, variant1, variant2, specific_mut_rates1,specific_mut_rates2) {

  tumors_with_variant1_mutated <- MAF_input1[MAF_input1$unique_variant_ID_AA==variant1,"Unique_patient_identifier"]
  tumors_with_variant2_mutated <- MAF_input2[MAF_input2$unique_variant_ID_AA==variant2,"Unique_patient_identifier"]

  tumors_with_both_mutated <- base::intersect(tumors_with_variant1_mutated,tumors_with_variant2_mutated)

  tumors_with_ONLY_variant1_mutated <- tumors_with_variant1_mutated[which(!tumors_with_variant1_mutated %in% tumors_with_variant2_mutated)]

  tumors_with_ONLY_variant2_mutated <- tumors_with_variant2_mutated[which(!tumors_with_variant2_mutated %in% tumors_with_variant1_mutated)]

  tumors_with_neither_mutated <- setdiff(rownames(all_tumors), c(tumors_with_both_mutated,tumors_with_ONLY_variant1_mutated,tumors_with_ONLY_variant2_mutated))

  if((gamma[3] == gamma[1] + gamma[2]) |
     (gamma[4] == gamma[1] + gamma[2])){return(-1e200)}

  sum_log_lik <- 0

  # not in either

  for (tumor in tumors_with_neither_mutated) {
    sum_log_lik <- sum_log_lik +
      (
        (-1*(gamma[1] * specific_mut_rates1[tumor, variant1])) +
          (-1*(gamma[2] * specific_mut_rates2[tumor, variant2]))
      )

  }


  # tumors only have variant 1


  for( tumor in tumors_with_ONLY_variant1_mutated){
    sum_log_lik <- sum_log_lik +

      # if(
      log(
        # log(
        -1*(
          (gamma[1] * specific_mut_rates1[tumor, variant1]) /
            (
              (gamma[1] * specific_mut_rates1[tumor, variant1]) +
                (gamma[2] * specific_mut_rates2[tumor, variant2]) -
                (gamma[4] * specific_mut_rates2[tumor, variant2])
            )
          # )
        ) *
          # log(
          ((exp(  (-1*(gamma[1] * specific_mut_rates1[tumor, variant1])) +
                    (-1*(gamma[2] * specific_mut_rates2[tumor, variant2])))) -
             (exp(-1*(gamma[4] * specific_mut_rates2[tumor, variant2]))))
        # )
      )
    # ==-Inf){break}

  }

  # tumors only have variant 2

  for( tumor in tumors_with_ONLY_variant2_mutated){
    sum_log_lik <- sum_log_lik +
      log(
        # log(
        -1* (
          (gamma[2] * specific_mut_rates2[tumor, variant2]) /
            (
              (gamma[2] * specific_mut_rates2[tumor, variant2]) +
                (gamma[1] * specific_mut_rates1[tumor, variant1]) -
                (gamma[3] * specific_mut_rates1[tumor, variant1])
            )
          # )
        ) *
          # log(
          ((exp(  (-1*(gamma[1] * specific_mut_rates1[tumor, variant1])) +
                    (-1*(gamma[2] * specific_mut_rates2[tumor, variant2])))) -
             (exp(-1*(gamma[3] * specific_mut_rates1[tumor, variant1]))))
        # )
      )

  }


  # tumors have both variants
  for( tumor in tumors_with_both_mutated) {

    sum_log_lik <- sum_log_lik + log(
      1-
        (
          ( # P(wt)
            exp((-1*(gamma[1] * specific_mut_rates1[tumor, variant1])) +
                  (-1*(gamma[2] * specific_mut_rates2[tumor, variant2])))
          ) +
            ( #P(1)
              # log(
              -1*(
                (gamma[1] * specific_mut_rates1[tumor, variant1]) /
                  (
                    (gamma[1] * specific_mut_rates1[tumor, variant1]) +
                      (gamma[2] * specific_mut_rates2[tumor, variant2]) -
                      (gamma[4] * specific_mut_rates2[tumor, variant2])
                  )
                # )
              ) *
                # log(
                ((exp(  (-1*(gamma[1] * specific_mut_rates1[tumor, variant1])) +
                          (-1*(gamma[2] * specific_mut_rates2[tumor, variant2])))) -
                   (exp(-1*(gamma[4] * specific_mut_rates2[tumor, variant2]))))
              # )
            ) +
            ( # P(2)
              -1* (
                (gamma[2] * specific_mut_rates2[tumor, variant2]) /
                  (
                    (gamma[2] * specific_mut_rates2[tumor, variant2]) +
                      (gamma[1] * specific_mut_rates1[tumor, variant1]) -
                      (gamma[3] * specific_mut_rates1[tumor, variant1])
                  )
                # )
              ) *
                # log(
                ((exp(  (-1*(gamma[1] * specific_mut_rates1[tumor, variant1])) +
                          (-1*(gamma[2] * specific_mut_rates2[tumor, variant2])))) -
                   (exp(-1*(gamma[3] * specific_mut_rates1[tumor, variant1]))))
              # )
            )


        )


    )



  }

  # in case it tried all the max at once.
  if(!is.finite(sum_log_lik)){
    return(-1e200)
  }
  return(sum_log_lik)
}
