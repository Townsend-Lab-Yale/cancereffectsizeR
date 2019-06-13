#' ml_objective
#'
#' Objective function that we will be optimizing in order to find the site specific selection intensity that maximizes the likelihood of each tumor having a mutation or not, where the mutation rates are site and tumor specific.
#'
#' @param gamma A selection intensity at which to calculate the likelihood
#' @param MAF_input A data frame that includes columns "Unique_patient_identifier", "Gene_name", and "unique_variant_ID_AA"
#' @param all_tumors A list of all the tumors we are calculating the likelihood across
#' @param stages An environment keys = tumor, value = all_tumors[tumor, "subset"]; or Null if all are in same class (subset = 1)
#' @param gene The gene we want to look at
#' @param variant The variant we want to look at
#' @param specific_mut_rates A matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants (produced by mutation_rate_calc)
#'
#' @return A log likelihood value
#' @export
#'
#' @examples





ml_objective <- function(gamma, MAF_input, all_tumors, stages, gene, variant, specific_mut_rates) {
  


  tumors_with_pos_mutated <- MAF_input$Unique_patient_identifier[MAF_input$unique_variant_ID_AA==variant & MAF_input$Gene_name==gene]
  tumors_without_gene_mutated <- rownames(all_tumors)[!rownames(all_tumors) %in% unique(MAF_input$Unique_patient_identifier[MAF_input$Gene_name==gene])]


  calc_gamma_sums = function(tumor) {
    if (is.null(stages)) {
      return(-1*gamma[1])
    }
    return(sum(-1*gamma[1:stages[[tumor]]]))
  }
  gamma_sums = vapply(tumors_without_gene_mutated, calc_gamma_sums, FUN.VALUE = 1.0)
  sum_log_lik = sum(gamma_sums * specific_mut_rates[tumors_without_gene_mutated, variant])


  ## this should be rewritten to match optimized code above for consistent style,
  ## but since any given position is mutated in very few tumors, the performance gains will minimal
  for (tumor in tumors_with_pos_mutated) {
    this_sum <- 0
    for(levels in 1:all_tumors[tumor,"subset"]){
      if(levels == 1){
        this_sum <- this_sum + 1-exp(-1*(gamma[levels] * specific_mut_rates[tumor, variant]))
      }else{
        level_of_occur <- 1-exp(-1*(gamma[levels] * specific_mut_rates[tumor, variant]))
        for(levels_sub in levels:2){
            level_of_occur <- level_of_occur * exp(-1*(gamma[levels_sub-1] * specific_mut_rates[tumor, variant]))
        }
        this_sum <- this_sum + level_of_occur
      }
    }
    sum_log_lik <- sum_log_lik + log(this_sum)
  }


  # in case it tried all the max at once.
  if(!is.finite(sum_log_lik)){
    return(-1e200)
  }
  return(sum_log_lik)
}
