#' cancereffectsizeR wrapper for fit_to_signatures 
#'
#' This function gets called internally by trinuc_mutation_rates() for each tumor in a CESAnalysis, accepting
#' a data.frame of mutation counts and returning fit_to_signatures output. Note: this function supports indels
#' if passed in the same format.
#' 
#' @param tumor_trinuc_counts matrix of trinuc variant counts where columns respond to tumors and 
#'   order of trinucleotide changes match signatures_df
#' @param signatures_df data.frame of signatures (see COSMIC v3 signatures included with package for format)
#' @param signatures_to_remove names of signatures in signatures_df to keep out of MutationalPatterns
#'   and assign zero weights. Only occurs when strict == FALSE
#' @param mp_strict_args named list of additional arguments to fit_to_signatures_strict
#' @param bootstrap_mutations T/F (default FALSE) whether to run
#'   fit_to_signatures_bootstrapped() with n_boot = 1 instead of fit_to_signatures_strict().
#' @keywords internal
#' @return a data.frame of signature weights
run_mutational_patterns = function(tumor_trinuc_counts, signatures_df, signatures_to_remove, 
                                   mp_strict_args = list(), bootstrap_mutations = FALSE) {

  signatures_to_include = signatures_df[! rownames(signatures_df) %in% signatures_to_remove,]
  
  # transpose so that columns correspond to signatures
  signatures_to_include = t(signatures_to_include)
  
  # convert signatures to matrix as required by MutationalPatterns
  signatures_to_include <- data.matrix(signatures_to_include)
  
  if(bootstrap_mutations) {
    mp_strict_args[['n_boot']] = 1
    if(! 'method' %in% mp_strict_args) {
      mp_strict_args[['method']] = 'strict'
    }
  }
  
  args = c(mp_strict_args, list(mut_matrix = tumor_trinuc_counts, signatures = signatures_to_include))
  
  # have to deal with a rare MP bug
  signatures_output = tryCatch(
  {
    if(bootstrap_mutations) {
      do.call(MutationalPatterns::fit_to_signatures_bootstrapped, args)
    } else {
      do.call(MutationalPatterns::fit_to_signatures_strict, args)$fit_res
    }
    
  }, error = function(e) {
    if (grepl('a dimension is zero', conditionMessage(e)))  {
      MutationalPatterns::fit_to_signatures(mut_matrix = tumor_trinuc_counts, signatures = signatures_to_include)
    } else {
      e
    }
  })
  
  # Convert bootstrap output to match fit_to_signatures_strict output,
  # and then produce a weights data.frame.
  if(bootstrap_mutations) {
    # In the bootstrap function, signatures with no weight get left out of output.
    # Add them to signatures_to_remove so they'll get put into table later.
    zeroed_sigs = setdiff(colnames(signatures_to_include), colnames(signatures_output))
    signatures_to_remove = union(zeroed_sigs, signatures_to_remove)
    weights = as.data.frame(signatures_output)
  } else {
    # must convert to data.frame and transpose so that result is compatible with
    # trinuc_mutation_rates
    weights = t(as.data.frame(signatures_output$contribution))
  }

  # add columns for any removed signatures into output matrix (with zero values)
  if(!is.null(signatures_to_remove)) {
    zeroed_sigs = as.data.frame(matrix(data = 0,nrow = 1,ncol = length(signatures_to_remove)))
    colnames(zeroed_sigs) <- signatures_to_remove
    weights <- cbind(weights, zeroed_sigs)
    
    # sort columns to match standard order
    weights <- weights[,rownames(signatures_df)]
  }
  
  return(weights)
}
