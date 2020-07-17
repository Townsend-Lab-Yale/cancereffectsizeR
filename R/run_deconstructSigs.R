#' cancereffectsizeR wrapper for deconstructSigs 
#'
#' This function gets called internally by trinuc_mutation_rates() for each tumor in a CESAnalysis, accepting its
#' a data.frame of mutation counts and returning deconstructSigs output.
#' 
#' Since this funciton implements artifact accounting and signature removal, you may find it useful for your own mutational 
#' signature analyses. However, you may be able to just use trinuc_mutation_rates().
#' 
#' @param tumor_trinuc_counts one-row data.frame of trinuc variant counts (in deconstructSigs order) for one tumor
#' @param signatures_df data.frame of signatures (see COSMIC v3 signatures included with package for format)
#' @param signatures_to_remove names of signatures in signatures_df to keep out of deconstructSigs and assign zero weights
#' @param artifact_signatures names of signatures that should be treated as experimental (e.g., sequencing) artifacts  
#'                            (these have weights calculated but are then normalized out when calculating true trinuc proportions)
#' @param tri.counts.method exome/genome trinucleotide content normalization argument to pass to deconstructSigs (see its docs)
#' @return a list with 2 items: raw_sig_output (exactly what comes out of deconstructSigs) and adjusted_sig_output, which contains
#'         signature weights (possibly adjusted with artifact accounting) and expected relative trinuc mutation rates based on 
#'         signature weights
#' @export
run_deconstructSigs = function(tumor_trinuc_counts, signatures_df, signatures_to_remove, 
                               artifact_signatures, tri.counts.method) {

  # toss all the "signatures_to_remove" from the complete set of signatures
  signatures_to_include = signatures_df[! rownames(signatures_df) %in% signatures_to_remove,]

  # contexts.needed indicates if trinuc normalization should happen; even if already normalized, 
  # harmless to normalize again
  signatures_output <- deconstructSigs::whichSignatures(tumor.ref = tumor_trinuc_counts,
                                                        signatures.ref = signatures_to_include,
                                                        contexts.needed = T,
                                                        tri.counts.method = tri.counts.method)


  # put raw dS output into output list
  output = list(raw_sig_output = signatures_output)

  
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

  
  # determine expected trinuc mutation proportions based on signature weights
  # in rare cases (possibly due to artifact accounting), all weights will be 0, and in that case just leave trinuc_prop NULL
  if(any(signatures_output$product != 0)) {
    trinuc_prop = signatures_output$product/sum(signatures_output$product)
    # Some trinuc SNVs have substitution rates of zero under certain signatures.
    # In rare cases, a tumor's fitted combination of signatures can therefore also
    # have a substitution rate of zero for particular trinucleotide contexts.
    # If this happens, we add the lowest nonzero rate to all rates and renormalize.
    if(any(trinuc_prop == 0)) {
      lowest_nonzero_rate = min(trinuc_prop[trinuc_prop != 0])
      trinuc_prop = trinuc_prop + lowest_nonzero_rate
      # renormalize so rates sum to 1
      trinuc_prop = trinuc_prop / sum(trinuc_prop)
    }
    signatures_output$trinuc_prop = trinuc_prop
  } else {
    signatures_output["trinuc_prop"] = list(NULL)
  }
  output$adjusted_sig_output = signatures_output[c("weights", "trinuc_prop")]
  return(output)
}