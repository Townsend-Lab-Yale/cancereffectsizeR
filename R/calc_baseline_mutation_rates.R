#' Runs functions supplying trinculeotide mutation weightings, mutation rate calculations, and MAF gene annotations
#' 
#'
#' @param cesa CESAnalysis object
#' @param covariate_file tissue-specific covariate file for dNdScv (gene-level mutation rate calculation)
#' @param trinuc_all_tumors Calculates trinucleotide signatures within all tumors (even those with < 50 variants)
#' @param signature_choice Either "signatures_cosmic_May2019" (default) or "signatures.cosmic" (COSMIC signatures v2 originally packaged with deconstructSigs).
#' @export



calc_baseline_mutation_rates <- function(
      cesa = NULL,
      covariate_file=NULL,
      #cores = 1, # currently unused, but could add multicore functionality for dNdScv and possibly deconstructSigs
      tumor_specific_rate_choice = F, # not currently used
      trinuc_all_tumors = T, # not currently used
      signature_choice = "signatures_cosmic_May2019",
      trinuc_algorithm_choice="weighted",
      ... ) { 

  # Calculate trinucleotide mutation weightings using deconstructSigs
  cesa = cancereffectsizeR::trinucleotide_mutation_weights(cesa, algorithm_choice = trinuc_algorithm_choice, signature_choice = signature_choice)

  # Calculate gene-level mutation rates using dNdScv
  cesa = cancereffectsizeR::gene_level_mutation_rates(cesa, covariate_file)

  # Assign genes to MAF, keeping assignments consistent with dndscv when possible
  cesa = cancereffectsizeR::annotate_gene_maf(cesa)


  return(cesa)
}


