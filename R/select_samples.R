#' Retrieve validated subset of CESAnalysis samples table
#' 
#' @param cesa CESAnalysis
#' @param samples Vector of Unique_Patient_Identifiers, or data.table consisting
#' of rows from a CESAnalysis samples table. If empty, returns full sample table.
#' @return data.table consisting of one or more rows from the CESAnalysis samples table.
#' @keywords internal
select_samples = function(cesa = NULL, samples = character()) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be CESAnalysis.")
  }
  if(! is(samples, "character") && ! is(samples, "data.table")) {
    stop("samples should be character or data.table.")
  }
  if (is(samples, "data.table")) {
    if(! "Unique_Patient_Identifier" %in% names(samples)) {
      stop("samples should be character or a data.table with a Unique_Patient_Identifier column.")
    }
    if(samples[, .N] == 0) {
      stop("The input samples table is empty")
    }
    samples = samples$Unique_Patient_Identifier
  }
  
  curr_sample_info = cesa@samples
  if(length(samples) > 0) {
    if(anyNA(samples)) {
      stop("NA values in samples.")
    }
    if(any(duplicated(samples))) {
      stop("Input samples contains duplicates.")
    }
    curr_sample_info = cesa@samples[samples, on = "Unique_Patient_Identifier", nomatch = NULL]
    if(curr_sample_info[, .N] != length(samples)) {
      stop("Some input samples not present in CESAnalysis samples table.")
    }
  }
  return(copy(curr_sample_info))
}