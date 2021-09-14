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


#' Find samples with specified variants
#'
#' A convenience function to identify samples with specific variants.
#'
#' @param cesa CESAnalysis
#' @param any_of Select samples with ANY of the given variant names/IDs, such as
#'   c("8:142506482_C>G", "KRAS G12C"). When a gene has multiple transcripts in reference
#'   data, you may wish to use full IDs, such as "KRAS_G12C_ENSP00000256078".
#' @export
samples_with = function(cesa, any_of = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa expected to be CESAnalysis.")
  }
  variants_by_type = sort_and_validate_variant_ids(cesa = cesa, input_ids = any_of, drop_unannotated = TRUE)
  
  snv_ids = variants_by_type[['snv_id']]
  aac_ids = variants_by_type[['aac_id']]
  
  snv_from_aac = unique(cesa@mutations$amino_acid_change[aac_ids, unlist(constituent_snvs), on = 'aac_id'])
  all_snv_ids = union(snv_from_aac, snv_ids)
  
  return(cesa@maf[all_snv_ids, unique(Unique_Patient_Identifier), on = 'variant_id', nomatch = NULL])
}
