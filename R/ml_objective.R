#' ml_objective
#'
#' Objective function that we will be optimizing in order to find the site specific selection intensity that maximizes the likelihood of each tumor having a mutation or not, where the mutation rates are site and tumor specific.
#'
#' @param gamma A selection intensity at which to calculate the likelihood
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
ml_objective <- function(gamma, MAF_input, all_tumors, gene, variant, specific_mut_rates) {

  tumors_with_pos_mutated <- MAF_input$Unique_patient_identifier[MAF_input$unique_variant_ID_AA==variant & MAF_input$Gene_name==gene]
  tumors_without_gene_mutated <- rownames(all_tumors)[!rownames(all_tumors) %in% unique(MAF_input$Unique_patient_identifier[MAF_input$Gene_name==gene])]



  sum_log_lik <- 0

  for (tumor in tumors_without_gene_mutated) {
    for(this_level in unique(all_tumors[,"subset"])) {
      sum_log_lik <- sum_log_lik + (-gamma[this_level] * specific_mut_rates[tumor, variant])
    }
  }



  for (tumor in tumors_with_pos_mutated) {
    for(levels in 1:all_tumors[tumor,"subset"]){
      sum_log_lik <- sum_log_lik + log(1-exp(-gamma[levels] * specific_mut_rates[tumor, variant]))
    }
  }



  #
  #   tumors_without_gene_mutated <- all_tumors[!all_tumors %in% unique(MAF_input$Unique_patient_identifier[MAF_input$Gene_name==gene])]
  #   tumors_with_pos_mutated <- MAF_input$Unique_patient_identifier[MAF_input$unique_variant_ID_AA==variant & MAF_input$Gene_name==gene]
  #
  #   sum_log_lik <- 0
  #   for (tumor in tumors_without_gene_mutated) {
  #     sum_log_lik <- sum_log_lik + (-gamma * specific_mut_rates[tumor, variant])
  #   }
  #   for (tumor in tumors_with_pos_mutated) {
  #     sum_log_lik <- sum_log_lik + log(1-exp(-gamma * specific_mut_rates[tumor, variant]))
  #   }


  return(sum_log_lik)
}
