#' Assign pre-calculated relative trinucleotide mutation rates
#'
#' This function assigns trinucleotide-context-specific relative SNV mutation rates to
#' tumors in a CESAnalysis. (These could be rates previously generated with
#' \code{trinuc_mutation_rates()}, or they could calculated using your own methods.) The
#' input rates must be a data.table or matrix. If supplying a data table, there must be a
#' Unique_Patient_Identifier column; if supplying a a matrix, the identifiers should be
#' supplied as rownames instead. Either way, all samples in the CESAnalysis must be
#' represented in the input rates. To avoid user error, there cannot be any superfluous
#' samples in the input rates unless \code{ignore_extra_samples = T}. Besides the
#' identifier column (or matrix rownames), there must be 96 columns, with column names
#' exactly matching the deconstructSigs/MutationalPatterns naming and order (run this
#' function with incorrect column names, and the names you need to use will be printed).
#' Since CES uses relative trinuc rates, rows must sum to 1, with all values greater than
#' 0. You'll get a warning if any rate is less than 1e-9, since (unrealistically) low
#' rates may crash selection model likelihood functions that aren't expecting such small
#' values.
#' 
#' 
#' @param cesa CESAnalysis object
#' @param trinuc_rates a matrix or data table (see description for format)
#' @param ignore_extra_samples skip samples in the input table that are not in the CESAnalysis (when false, will stop with an error)
#' @export
set_trinuc_rates = function(cesa, trinuc_rates, ignore_extra_samples = FALSE) {
  if(! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be a CESAnalysis object", call. = F)
  }
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  if(cesa@samples[, .N] == 0) {
    stop("There are no samples in the CESAnalysis", call. = F)
  }
  if(is.null(trinuc_rates) || ! is(trinuc_rates, "matrix")) {
    if(is(trinuc_rates, "data.table")) {
      if(! "Unique_Patient_Identifier" %in% names(trinuc_rates)) {
        stop("If you supply a data.table, there must be a Unique_Patient_Identifier column", call. = F)
      }
      trinuc_proportion_matrix = as.matrix(trinuc_rates[, -"Unique_Patient_Identifier"])
      rownames(trinuc_proportion_matrix) = trinuc_rates$Unique_Patient_Identifier
    } else {
      stop("Expected trinuc rates to be a data.table or matrix", call. = F)
    }
  } else {
    trinuc_proportion_matrix = trinuc_rates
  }
  
  if(anyNA(trinuc_proportion_matrix)) {
    stop("NA values found in trinuc_proportion_matrix", call. = F)
  }
  
  if(ncol(trinuc_proportion_matrix) != 96 || ! identical(colnames(trinuc_proportion_matrix), deconstructSigs_trinuc_string)) {
    message("Expected column names:")
    names_with_quotes = paste0('"', deconstructSigs_trinuc_string, '"')
    cat(names_with_quotes, sep = ',')
    cat("\n")
    stop("Incorrect number of columns and/or incorrect column names or column order", call. = F)
  }
  
  if(length(setdiff(rownames(trinuc_proportion_matrix), cesa$samples$Unique_Patient_Identifier)) > 0) {
    if(! ignore_extra_samples) {
      stop(paste0("There are samples in the input rates that are not in the CESAnalysis; if this is intentional, re-run with \n",
                  "ignore_extra_samples = TRUE."), call. = F)
    }
    trinuc_proportion_matrix = trinuc_proportion_matrix[cesa$samples$Unique_Patient_Identifier, , drop = F]
  }
  if(! identical(sort(cesa@samples$Unique_Patient_Identifier), sort(rownames(trinuc_proportion_matrix)))) {
    stop("Not all samples in the CESAnalysis are present in the input rates.", call. = F)
  }
  
  
  # rows must sum to 1, but allow small tolerance 
  if(any(abs(rowSums(trinuc_proportion_matrix) - 1) > .001)) {
    stop("row sums of input rates must all be 1", call. = F)
  }
  if (any(apply(trinuc_proportion_matrix, 1, function(x) any(x <= 1e-9)))) {
    # rate < 0 is invalid
    if (any(apply(trinuc_proportion_matrix, 1, function(x) any(x < 0)))) {
      stop("input rates cannot have negative entries.")
    }
    msg = paste0("Some relative mutation rates are very low (<1e-9). This could cause problems, especially if",
                   " any \"impossible\" mutations (rate = 0) are in the data set. Very low rates, besides being unrealistic, can also ",
                   "crash selection model likelihood functions due to numerical precision issues. Consider tweaking your method!")
    warning(pretty_message(msg, emit = F))
  }
  cesa@trinucleotide_mutation_weights = list(trinuc_proportion_matrix = trinuc_proportion_matrix)
  cesa@advanced$locked = T
  cesa@advanced$trinuc_done = T
  return(cesa)
}
    
    