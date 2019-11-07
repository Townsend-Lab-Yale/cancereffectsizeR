#' Runs functions supplying trinculeotide mutation weightings, mutation rate calculations, and MAF gene annotations
#'
#'
#' @param cesa CESAnalysis object
#' @param covariate_file tissue-specific covariate file for dNdScv (gene-level mutation rate calculation)
#' @param signature_choice Either "signatures_cosmic_May2019" (default) or "signatures.cosmic" (COSMIC signatures v2 originally packaged with deconstructSigs).
#' @param signatures_to_remove Removes signatures from full potential signatures to minimize "signature bleeding", along the rationale proposed within the manuscript that originally calculated the v3 signature set doi: https://doi.org/10.1101/322859. Use `NULL` to keep all signatures. Signatures must correspond to the signature names in `signature_choice`.
#' @param trinuc_algorithm_choice
#' @param artifact_accounting
#' @param BIG_mem if FALSE then it deletes some information that, while maybe useful to the user, is unnecessary for the `effect_size_SNV` function
#' @param ...
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
      mutation_count_rules = T,
      BIG_mem = TRUE,
      ... ) {

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

  if(!BIG_mem){
    levels_in_selection_analysis <- names(cesa@progressions@order)

      for(subset_index in 1:length(levels_in_selection_analysis)){

       cesa@dndscv_out_list[[subset_index]] <-  list(sel_cv = cesa@dndscv_out_list[[subset_index]]$sel_cv)
      }


      # only need RefCDS that will be useful downstream

      # cesa@refcds_data <- cesa@refcds_data[unique(cesa@annotated.snv.maf$Gene_name[which(cesa@annotated.snv.maf$next_to_splice == T)])]


      list_extract <- function(x){
        return(list(gene_name=x$gene_name,
                    gene_id = x$gene_id,
                    seq_cds=x$seq_cds,
                    seq_cds1up=x$seq_cds1up,
                    seq_cds1down=x$seq_cds1down))
      }

      cesa@refcds_data <- as.array(lapply(cesa@refcds_data, list_extract))


      # list_extract <- function(x){
      #   return(list(gene_name=x$gene_name,
      #               gene_id = x$gene_id))
      # }


      genes_to_keep_info <- unique(cesa@annotated.snv.maf$Gene_name[which(cesa@annotated.snv.maf$next_to_splice == T)])

      for(gene_ind in 1:length(cesa@refcds_data)){
        if(!cesa@refcds_data[[gene_ind]]$gene_name %in% genes_to_keep_info){
          cesa@refcds_data[[gene_ind]] <- as.array(list(gene_name=cesa@refcds_data[[gene_ind]]$gene_name,
                                                        gene_id = cesa@refcds_data[[gene_ind]]$gene_id))
        }
      }




    }
  return(cesa)
}


