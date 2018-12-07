# Objective function that we will be optimizing in order to find the site specific selection
# intensity that maximizes the likelihood of each tumor having a mutation or not, where the
# mutation rates are site and tumor specific.

# Inputs:
# gamma: a selection intensity at which to calculate the likelihood
# MAF_input: a data frame that includes columns "Unique_patient_identifier", "Gene_name", and "unique_variant_ID_AA"
# all_tumors: a list of all the tumors we are calculating the likelihood across
# gene, variant: the gene and variant we want to look at
# specific_mut_rates: a matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants

# Outputs:
# A log likelihood value


ml_objective <- function(gamma, MAF_input, all_tumors, gene, variant, specific_mut_rates) {
  tumors_without_gene_mutated <- all_tumors[!all_tumors %in% unique(MAF_input$Unique_patient_identifier[MAF_input$Gene_name==gene])]
  tumors_with_pos_mutated <- MAF_input$Unique_patient_identifier[MAF_input$unique_variant_ID_AA==variant & MAF_input$Gene_name==gene]

  sum_log_lik <- 0
  for (tumor in tumors_without_gene_mutated) {
    sum_log_lik <- sum_log_lik + (-gamma * specific_mut_rates[tumor, variant])
  }
  for (tumor in tumors_with_pos_mutated) {
    sum_log_lik <- sum_log_lik + log(1-exp(-gamma * specific_mut_rates[tumor, variant]))
  }

  return(sum_log_lik)
}
