#' Runs functions supplying trinculeotide mutation weightings, mutation rate calculations, and MAF gene annotations
#'
#'
#' @param cesa CESAnalysis object
#' @param covariate_file tissue-specific covariate file for dNdScv (gene-level mutation rate calculation)
#' @param signature_choice Either "signatures_cosmic_May2019" (default) or "signatures.cosmic" (COSMIC signatures v2 originally packaged with deconstructSigs).
#' @param signatures_to_remove Removes signatures from full potential signatures to minimize "signature bleeding", along the rationale proposed within the manuscript that originally calculated the v3 signature set doi: https://doi.org/10.1101/322859. Use `NULL` to keep all signatures. Signatures must correspond to the signature names in `signature_choice`.
#' @param trinuc_algorithm_choice
#' @param artifact_accounting
#' @param mutation_count_rules T/F on whether to follow the mutation count rules outlined in https://doi.org/10.1101/322859, the manuscript reported the v3 COSMIC signature set.
#'
#' @export



calc_baseline_mutation_rates <- function(
      cesa = NULL,
      covariate_file=NULL,
      #cores = 1, # currently unused, but could add multicore functionality for dNdScv and possibly deconstructSigs
      signature_choice = "signatures_cosmic_May2019",
      trinuc_algorithm_choice="weighted",
      artifact_accounting = T,
      signatures_to_remove = c("SBS25","SBS31","SBS32","SBS35"),
      mutation_count_rules = T) 
{

  # Calculate trinucleotide mutation weightings using deconstructSigs
  cesa = cancereffectsizeR::trinucleotide_mutation_weights(cesa,
                                                           algorithm_choice = trinuc_algorithm_choice,
                                                           signature_choice = signature_choice,
                                                           artifact_accounting = artifact_accounting,
                                                           signatures_to_remove = signatures_to_remove)

  # Calculate gene-level mutation rates using dNdScv
  cesa = cancereffectsizeR::gene_level_mutation_rates(cesa, covariate_file)

  # Assign genes to MAF, keeping assignments consistent with dndscv when possible
  cesa = cancereffectsizeR::annotate_gene_maf(cesa)


  return(cesa)
}


