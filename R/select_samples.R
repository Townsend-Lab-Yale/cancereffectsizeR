#' Retrieve validated subset of CESAnalysis samples table
#' 
#' @param cesa CESAnalysis
#' @param samples Vector of patient_ids, or data.table consisting
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
    if(! "patient_id" %in% names(samples)) {
      stop("samples should be character or a data.table with a patient_id column.")
    }
    if(samples[, .N] == 0) {
      stop("The input samples table is empty")
    }
    samples = samples$patient_id
  }
  
  curr_sample_info = cesa@samples
  if(length(samples) > 0) {
    if(anyNA(samples)) {
      stop("NA values in samples.")
    }
    if(any(duplicated(samples))) {
      stop("Input samples contains duplicates.")
    }
    curr_sample_info = cesa@samples[samples, on = "patient_id", nomatch = NULL]
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
  
  if(is(any_of, 'VariantSetList')) {
    stop('any_of should be type character. For a VariantSetList, use [VariantSetList]$samples_with.')
  }
  variants_by_type = sort_and_validate_variant_ids(cesa = cesa, input_ids = any_of, drop_unannotated = TRUE)
  
  sbs_ids = variants_by_type[['sbs_id']]
  aac_ids = variants_by_type[['aac_id']]
  
  sbs_from_aac = cesa@mutations$aac_sbs_key[aac_ids, unique(sbs_id), on = 'aac_id']
  all_sbs_ids = union(sbs_from_aac, sbs_ids)
  
  return(cesa@maf[all_sbs_ids, unique(patient_id), on = 'variant_id', nomatch = NULL])
}
