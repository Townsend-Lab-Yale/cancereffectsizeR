#' deconstructSigs wrapper
#'
#' This function gets called internally during tumor-specific trinucleotide mutation rate calculation.
#' 
#' @param tumor_trinuc_counts one-row data.frame of trinuc variant counts (in deconstructSigs order) for one tumor
#' @param signatures_df data.frame of signatures (see COSMIC v3 signatures included with package for format)
#' @param signatures_to_remove names of signatures in signatures_df to keep out of deconstructSigs and assign zero weights
#' @param artifact_signatures names of signatures that should be treated as experimental (e.g., sequencing) artifacts  
#'                            (these have weights calculated but are then normalized out when calculating true trinuc proportions)
#' @param tri.counts.method exome/genome trinucleotide content normalization argument to pass to deconstructSigs (see its docs)
#' @export
run_deconstructSigs = function(tumor_trinuc_counts, signatures_df, signatures_to_remove, 
                               artifact_signatures, tri.counts.method) {

  # toss all the "signatures_to_remove" from the complete set of signatures
  signatures_to_include = signatures_df[! rownames(signatures_df) %in% signatures_to_remove,]

  # contexts.needed indicates if trinuc normalization should happen
  contexts.needed = ifelse(is(tri.counts.method, "character") && tri.counts.method == "default", FALSE, TRUE)
  
  signatures_output <- deconstructSigs::whichSignatures(tumor.ref = tumor_trinuc_counts,
                                                        signatures.ref = signatures_to_include,
                                                        contexts.needed = contexts.needed,
                                                        tri.counts.method = tri.counts.method)


  # add columns for any removed signatures into output matrix (with zero values)
  if(!is.null(signatures_to_remove)){
    zeroed_sigs = as.data.frame(matrix(data = 0,nrow = 1,ncol = length(signatures_to_remove)))
    colnames(zeroed_sigs) <- signatures_to_remove
    signatures_output$weights <- cbind(signatures_output$weights, zeroed_sigs)
    
    # sort columns to match standard order
    signatures_output$weights <- signatures_output$weights[,rownames(signatures_df)]
  }

  # We remove artifact signatures and renormalize so that we can
  # determine mutational flux in tumor from true sources
  if(! is.null(artifact_signatures)) {
    if(any(signatures_output$weights[artifact_signatures] > 0 )){
      # often, sum of best estimation of signatures is < 1.
      # will renormalize to this original total explanation of weight.
      initial_weight_sum <- sum(signatures_output$weights)

      # remove artifacts from mutational signatures flux
      signatures_output$weights[artifact_signatures] <- 0

      # if other weights exist, renormalize; otherwise weights get zeroed out (unlikely)
      if(any(signatures_output$weights>0)){
        signatures_output$weights <- signatures_output$weights/(sum(signatures_output$weights)/initial_weight_sum)
      }else{
        signatures_output$weights <- signatures_output$weights * 0
      }
      signatures_output$product <- as.matrix(signatures_output$weights) %*% as.matrix(signatures_df)
    }
  }
  
  # return signatures_output, a list which should look like normal deconstructSigs output
  return(signatures_output)
}