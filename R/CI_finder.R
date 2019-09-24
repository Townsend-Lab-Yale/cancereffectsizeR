#' MLE CI finder
#'
#' @param gamma_max selection intensity that maximizes the objective function
#' @param MAF_input A data frame that includes columns "Unique_Patient_Identifier", "Gene_name", and "unique_variant_ID_AA"
#' @param all_tumors A list of all the tumors we are calculating the likelihood across
#' @param progressions CESProgressions object
#' @param gene The gene we want to look at
#' @param variant The variant we want to look at
#' @param specific_mut_rates A matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants (produced by mutation_rate_calc)
#' @param log_units_down log units below maximum that correspond to 99.9 percent confidence
#' @param threshold tolerance threshold for the uniroot function
#'
#' @return
#' @export
#'
#' @examples
CI_finder <- function(gamma_max,
                      all_tumors,
                      MAF_input,
                      progressions,
                      gene,
                      variant,
                      specific_mut_rates,
                      log_units_down = 10.828/2,
                      threshold = 1e-2){



  # Temporary fix (simplify CI_finder parameters later)
  tumors_with_pos_mutated <- MAF_input$Unique_Patient_Identifier[MAF_input$unique_variant_ID_AA==variant & MAF_input$Gene_name==gene]
  tumors_without_gene_mutated <- all_tumors[ ! all_tumors %in% unique(MAF_input$Unique_Patient_Identifier[MAF_input$Gene_name==gene])]
  tumor_stages = get_progression_number(progressions, all_tumors)


  max_likelihood <- cancereffectsizeR::ml_objective(gamma = gamma_max, tumor_stages = tumor_stages, tumors_without_gene_mutated = tumors_without_gene_mutated,
      tumors_with_pos_mutated = tumors_with_pos_mutated, variant=variant, specific_mut_rates=specific_mut_rates)

  max_likelihood_CI <- max_likelihood - log_units_down



  #find lower CI

  lower_answer <- uniroot(f = cancereffectsizeR::ml_objective , lower = 0,upper = gamma_max, tumor_stages = tumor_stages, tumors_without_gene_mutated = tumors_without_gene_mutated,
      tumors_with_pos_mutated = tumors_with_pos_mutated, variant=variant, specific_mut_rates=specific_mut_rates, modifier = max_likelihood_CI,tol = threshold)


  #find upper CI

  upper_answer <- uniroot(f = cancereffectsizeR::ml_objective , lower = gamma_max,upper = 1e20, tumor_stages = tumor_stages, tumors_without_gene_mutated = tumors_without_gene_mutated,
      tumors_with_pos_mutated = tumors_with_pos_mutated, variant=variant, specific_mut_rates=specific_mut_rates, modifier = max_likelihood_CI,tol = threshold)


  return(list(lower_CI=lower_answer$root,
              upper_CI=upper_answer$root))

}
