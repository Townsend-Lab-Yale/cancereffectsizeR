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
  
  if(! is.character(any_of)) {
    stop("any_of should be character vector of variant identifiers")
  }
  any_of = unique(na.omit(any_of))
  if(length(any_of) == 0) {
    stop("any_of doesn't list any variants.")
  }
  snv_ids = intersect(any_of, cesa@mutations$snv$snv_id)
  any_of = setdiff(any_of, snv_ids)
  
  aac_ids = intersect(any_of, cesa@mutations$amino_acid_change$aac_id)
  any_of = setdiff(any_of, aac_ids)
  
  ## To-do: insert DBS/indel logic here (and below)
  
  if(length(any_of) > 0) {
    any_of = sub(' ', '_', any_of)
    tmp_dt = copy(cesa@mutations$amino_acid_change)
    tmp_dt[, tmp_name :=  paste(gene, aachange, sep = '_')]
    more_aac_ids = tmp_dt[any_of, aac_id, on = 'tmp_name', nomatch = NULL]
    any_of = setdiff(any_of, tmp_dt[more_aac_ids, tmp_name, on = 'aac_id'])
    aac_ids = union(aac_ids, more_aac_ids)
  }
  
  # Keep going if all unmatched IDs look valid. (Presumably, these are absent from annotations.)
  if(length(any_of) > 0) {
    apparent_snvs = any_of[grepl(':\\d+_[ACTG]>[ACGT]$', any_of)]
    if(length(apparent_snvs) > 0) {
      validate_snv_ids(apparent_snvs, get_cesa_bsg(cesa))
    }
    # insert more indel/dbs logic
    apparent_aac = setdiff(any_of, apparent_snvs)
    apparent_aac = complete_aac_ids(apparent_aac, .ces_ref_data[[cesa@ref_key]])
    aac_problems = validate_aac_ids(apparent_aac, .ces_ref_data[[cesa@ref_key]])
    if(! is.null(aac_problems)) {
      print(aac_problems)
      stop("The above variant IDs from any_of appear invalid.")
    }
  }
  snv_from_aac = unique(cesa@mutations$amino_acid_change[aac_ids, unlist(constituent_snvs), on = 'aac_id'])
  all_snv_ids = union(snv_from_aac, snv_ids)
  
  return(cesa@maf[all_snv_ids, unique(Unique_Patient_Identifier), on = 'variant_id', nomatch = NULL])
}

