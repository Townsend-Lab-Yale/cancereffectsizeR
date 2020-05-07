#' Assign pre-calculated relative trinucleotide mutation rates
#'
#' This function assigns a matrix of trinucleotide-context-specific relative 
#' SNV mutation rates to a CESAnalysis. You can use this function to restore 
#' previously calculated rates or to provide rates calculated using your own 
#' methods (i.e., not trinuc_mutation_rates() and assign_group_average_trinuc_rates()).
#' The maxtrix must cover all samples in the CESAnalysis, with row names matching
#' sample identifiers. There must be 96 columns, with column names exactly matching
#' the deconstructSigs naming and order (run this function with incorrect column names, 
#' and the names you need to use will be printed). Since CES uses relative trinuc rates,
#' rows must sum to 1, with all values greater than 0.
#' 
#' @param cesa CESAnalysis object
#' @param trinuc_proportion_matrix rows=tumors, columns=relative mutation rates (see description)
#' @export
set_trinuc_rates = function(cesa, trinuc_proportion_matrix) {
  if(is.null(cesa) || ! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be a CESAnalysis object")
  }
  if(cesa@maf[, .N] == 0) {
    stop("No MAF data in the CESAnalysis")
  }
  if(is.null(trinuc_proportion_matrix) || ! is(trinuc_proportion_matrix, "matrix")) {
    stop("Expected trinuc_proportion_matrix to be a matrix")
  }
  
  if(anyNA(trinuc_proportion_matrix)) {
    stop("NA values found in trinuc_proportion_matrix")
  }
  
  if(ncol(trinuc_proportion_matrix) != 96 || ! identical(colnames(trinuc_proportion_matrix), deconstructSigs_trinuc_string)) {
    message("Expected column names:")
    names_with_quotes = paste0('"', deconstructSigs_trinuc_string, '"')
    cat(names_with_quotes, sep = ',')
    cat("\n")
    stop("Incorrect number of columns and/or incorrect column names or column order")
  }
  
  if(! identical(sort(cesa@samples$Unique_Patient_Identifier), sort(rownames(trinuc_proportion_matrix)))) {
    stop("Row names of trinuc_proprotion_matrix must exactly match sample names (with none extra or missing); see samples([CESAnalysis])")
  }
  
  
  # rows must sume to 1, but allow small tolerance 
  if(any(abs(rowSums(trinuc_proportion_matrix) - 1) > .001)) {
    stop("row sums of matrix must all be 1")
  }
  if (any(apply(trinuc_proportion_matrix, 1, function(x) any(x <= 0)))) {
    # rate < 0 is invalid
    if (any(apply(trinuc_proportion_matrix, 1, function(x) any(x < 0)))) {
      stop("matrix cannot have negative entries.")
    }
    warning(paste0("Some mutation rates are zero in some tumors. This could cause problems, especially if",
                   "any \"impossible\" mutations are in the data set. Consider tweaking your method!"))
  }
  cesa@trinucleotide_mutation_weights = list(trinuc_proportion_matrix = trinuc_proportion_matrix)
  cesa@status[["trinucleotide mutation rates"]] = "User-supplied via set_trinuc_rates"
  cesa@status[["gene mutation rates"]] = "uncalculated (run gene_mutation_rates)"
  return(cesa)
}
    
    