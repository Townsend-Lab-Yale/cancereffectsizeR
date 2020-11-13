#' Set SNV signature weights
#' 
#' If you wish to use your own method to calculate sample-specific SNV signature weights
#' (as opposed the signature extraction built into trinuc_mutation_rates()), you can use
#' this function to load them into the CESAnalysis. Your input signatures will be used to
#' infer relative trinucleotide-context-specific mutation rates for all tumors. (This
#' means you can run set_signature_weights() or set_trinuc_rates(), but not both.) As in
#' trinuc_mutation_rates(), you can use a built-in set of signatures, such as COSMIC_v3.1,
#' or you can supply your own signature set definitions as documented in
#' ?trinuc_mutation_rates.
#'
#' The input data table must have a Unique_Patient_Identifier column and one column per
#' signature in the signature set. All samples in the CESAnalysis must be included in the
#' input table, and each sample's weights should have a sum on (0, 1]. Since these weights
#' are used by cancerefffectsizeR to infer trinucleotide-context-specific relative rates
#' of SNV mutations, each sample must have at least one non-artifact signature with
#' nonzero weight. (In the unlikely event that this is a problem, consider assigning
#' group-average signature weights to the artifact-only samples.)
#' 
#' @param cesa CESAnalysisy
#' @param signature_set signature set name (see \code{list_ces_signature_sets()}), or a
#'   custom signature set (see documentation in \code{trinuc_mutation_rates()})
#' @param weights data.table of relative signature weights for each sample (see details)
#' @param ignore_extra_samples skip samples in the input table that are not in the
#'   CESAnalysis (when false, will stop with an error)
#' @export
set_signature_weights = function(cesa, signature_set, weights, ignore_extra_samples = FALSE) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be CESAnalysi")
  }
  
  if (length(cesa@trinucleotide_mutation_weights) > 0) {
    stop("Trinucleotide mutation rates have already been calculated for this CESAnalysis, so you can't use set_signature_weights.")
  }
  
  if(is(signature_set, "character")) {
    if (length(signature_set) != 1) {
      stop("signature_set should be 1-length character; run list_ces_signature_sets() for options.\n(Or, it can a custom signature set; see docs.)", call. = F)
    }
    signature_set_data = get_ces_signature_set(cesa@ref_key, signature_set)
  } else if(is(signature_set, "list")) {
    validate_signature_set(signature_set)
    signature_set_data = signature_set
    check = tryCatch(get_ces_signature_set(cesa@ref_key,  signature_set_data$name),
                     error = function(e) NULL)
    if (! is.null(check)) {
      stop("Your signature set's name matches one already provided with your CESAnalysis reference data.\n",
           "If you're trying to supply a custom signature set, change its name. Otherwise, just specify the set by name.")
    }
  } else {
    stop("signature_set should be type character; run list_ces_signature_sets() for options.\n",
         "(Or, it can be a custom signature set; see details in ?trinuc_mutation_rates().)")
  }
  signature_set_name = signature_set_data$name
  signatures = signature_set_data$signatures
  signature_metadata = signature_set_data$meta
  
  if(! is(weights, "data.table")) {
    stop("weights should be a data.table (you may need to convert with as.data.table()).")
  }
  if (! is.logical(ignore_extra_samples) || length(ignore_extra_samples) != 1) {
    stop("ignore_extra samples should be T/F")
  }
  
  if (uniqueN(colnames(weights)) != length(weights)) {
    stop("Input weights table has repeated column names.")
  }
  
  # validate weights
  signature_names = rownames(signatures)
  if (! identical(sort(c(signature_names, "Unique_Patient_Identifier")), sort(colnames(weights)))) {
    stop("Column names of weights must exactly match all signature names in the signature set (with none missing), ",
         "plus a Unique_Patient_Identifier column.")
  }
  
  if (! all(weights[, sapply(.SD, is.numeric), .SDcols = signature_names])) {
    stop("Not all signature weight columns in input are type numeric.")
  }
  if (! is(weights$Unique_Patient_Identifier, "character")) {
    stop("Input weights column Unique_Patient_Identifier should be type character.")
  }
  
  if(any(duplicated(weights$Unique_Patient_Identifier))) {
    stop("Some samples appear more than once in the weights table.")
  }
  
  
  if (length(setdiff(cesa@samples$Unique_Patient_Identifier, weights$Unique_Patient_Identifier)) > 0) {
    stop("Not all samples in the CESAnalysis appear in the input weights table.")
  }
  
  if (length(setdiff(weights$Unique_Patient_Identifier, cesa@samples$Unique_Patient_Identifier)) > 0 & ! ignore_extra_samples) {
    stop("There are samples in the input weight table that are not part of the CESAnalysis. (Re-run with ",
         "ignore_extra_samples = TRUE to ignore these.)")
  }
  
  # subset to CESAnalysis samples
  weights = weights[cesa@samples$Unique_Patient_Identifier, on = "Unique_Patient_Identifier"]
  
  
  if(anyNA(weights)) {
    stop("There are NA values in the input weights table.")
  }
  
  weights_okay = apply(weights[, -"Unique_Patient_Identifier"], 1, function(x) { total = sum(x); total > 0 && total <= 1})
  if (! all(weights_okay)) {
    bad_samples = weights[which(! weights_okay), Unique_Patient_Identifier]
    stop("Some samples have all-zero weights or weights summing greater than 1:\n", paste(bad_samples, collapse = ", "), ".")
  }

  
  total_snv_counts = cesa@maf[variant_type == "snv"][, .(total_snvs = .N), by = "Unique_Patient_Identifier"]
  weights[total_snv_counts, total_snvs := total_snvs, on = "Unique_Patient_Identifier"]
  weights[is.na(total_snvs), total_snvs := 0]
  weights[, c("sig_extraction_snvs", "group_avg_blended") := list(NA_integer_, NA)]
  setcolorder(weights, c("Unique_Patient_Identifier", "total_snvs", "sig_extraction_snvs", "group_avg_blended"))
  
  # artifact accounting
  bio_weights = copy(weights)
  if (sum(signature_metadata$Likely_Artifact) > 0) {
    artifact_signatures = signature_metadata[Likely_Artifact == TRUE, Signature]
    bio_weights[, initial_weight_sum := rowSums(.SD), .SDcols = signature_names]
    bio_weights[, (artifact_signatures) := 0]
    zeroed_out_index = which(bio_weights[, rowSums(.SD), .SDcols = signature_names] == 0)
    if (length(zeroed_out_index) > 0) {
      stop("Some tumor(s) have all-zero weights across non-artifact signatures:\n",
              paste(bio_weights[zeroed_out_index, Unique_Patient_Identifier], collapse = ", "), ".")
    }
    bio_weights[zeroed_out_index, adjust := 1]
    bio_weights[! zeroed_out_index, adjust := sum(.SD)/initial_weight_sum, by = "Unique_Patient_Identifier", .SDcols = signature_names]
    bio_weights[, (signature_names) := lapply(.SD, `/`, adjust), .SDcols = signature_names]
    bio_weights[, c("initial_weight_sum", "adjust") := NULL]
  }
  
  # get trinuc rates and normalize to relative rates
  trinuc_rates = as.matrix(bio_weights[, ..signature_names]) %*% as.matrix(signatures)
  trinuc_rates = t(apply(trinuc_rates, 1, function(x) x / sum(x))) # apply transposes against our wishes
  rownames(trinuc_rates) = bio_weights$Unique_Patient_Identifier
  
  # if any sample has a rate of zero in any context (possible with the right combination of signatures),
  # add the lowest nonzero rate to all and renormalze
  add_pseudo = function(x) {
    if(any(x == 0)) {
      next_lowest = min(x[x != 0])
      x = x + next_lowest
      x = x / sum(x)
    }
    return(x)
  }
  trinuc_rates = t(apply(trinuc_rates, 1, add_pseudo))
  
  cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix = trinuc_rates
  cesa@trinucleotide_mutation_weights$signature_weight_table_with_artifacts = weights
  cesa@trinucleotide_mutation_weights$signature_weight_table = bio_weights
  cesa@advanced$locked = TRUE
  cesa@advanced$trinuc_done = TRUE
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}