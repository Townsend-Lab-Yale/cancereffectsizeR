#' cancereffectsizeR wrapper for deconstructSigs 
#'
#' This function gets called internally by trinuc_mutation_rates() for each tumor in a CESAnalysis, accepting
#' a data.frame of mutation counts and returning a data.frame of signature weights.
#' 
#' @param tumor_trinuc_counts one-row data.frame of trinuc variant counts (in deconstructSigs order) for one tumor
#' @param signatures_df data.frame of signatures (see COSMIC v3 signatures included with package for format)
#' @param signatures_to_remove names of signatures in signatures_df to keep out of deconstructSigs and assign zero weights
#' @param tri.counts.method exome/genome trinucleotide content normalization argument to pass to deconstructSigs (see its docs)
#' @keywords internal
#' @return a data.frame of signature weights
run_deconstructSigs = function(tumor_trinuc_counts, signatures_df, signatures_to_remove, tri.counts.method) {
  # toss all the "signatures_to_remove" from the complete set of signatures
  signatures_to_include = signatures_df[! rownames(signatures_df) %in% signatures_to_remove,]

  # contexts.needed indicates if trinuc normalization should happen; even if already normalized, 
  # harmless to normalize again
  signatures_output <- deconstructSigs::whichSignatures(tumor.ref = tumor_trinuc_counts,
                                                        signatures.ref = signatures_to_include,
                                                        contexts.needed = T,
                                                        tri.counts.method = tri.counts.method)

  all_weights = signatures_output$weights
  
  if (!is.null(signatures_to_remove)) {
    zeroed_sigs = as.data.frame(matrix(data = 0,nrow = 1,ncol = length(signatures_to_remove)))
    colnames(zeroed_sigs) <- signatures_to_remove
    all_weights <- cbind(all_weights, zeroed_sigs)
    
    # sort columns to match standard order
    all_weights = all_weights[,rownames(signatures_df)]
  }
  
  return(all_weights)
}