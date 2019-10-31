#' Calculate selection intensity
#'
#' Actual function to find the site specific selection intensity that maximizes the likelihood of each tumor being mutated or not. Uses site and tumor specific mutation rates. Uses Brent 1 dimensional optimization technique.optimize_gamma
#'
#' @param MAF_input A data frame that includes columns "Unique_patient_identifier", "Gene_name", and "unique_variant_ID_AA"
#' @param eligible_tumors a vector of tumors we are calculating the likelihood across (excludes tumors without coverage at the variant site)
#' @param gene The gene we want to look at
#' @param progressions CESProgressions object
#' @param variant The variant we want to look at
#' @param specific_mut_rates A matrix of site and tumor specific mutation rates where the rows correspond to tumors and the columns to variants (produced by mutation_rate_calc)
#'
#' @return The optimal selection intensity for the gene and variant that maximizes the likelihood of the observations
#' @export
#'
#' @examples

optimize_gamma <- function(MAF_input, eligible_tumors, progressions, gene, variant, specific_mut_rates) {
  par_init <- rep(1000,length=length(progressions@order))
  tumors_with_pos_mutated <- MAF_input$Unique_Patient_Identifier[MAF_input$unique_variant_ID_AA==variant & MAF_input$Gene_name==gene]
  # specific_mut_rates contains only tumors with some mutations in data set (some tumors with coverage may not have any mutations), so subset
  tumors_without_gene_mutated <- eligible_tumors[eligible_tumors %in% rownames(specific_mut_rates) & ! eligible_tumors %in% unique(MAF_input$Unique_Patient_Identifier[MAF_input$Gene_name==gene])]
  tumor_stages = get_progression_number(progressions, eligible_tumors)
    return(optim(par=par_init, fn=ml_objective, tumor_stages = tumor_stages, tumors_without_gene_mutated = tumors_without_gene_mutated,
    	tumors_with_pos_mutated = tumors_with_pos_mutated, variant=variant, specific_mut_rates=specific_mut_rates,
                 method="L-BFGS-B", lower=1e-3, upper=1000000000, control=list(fnscale=-1e-12)))
}
