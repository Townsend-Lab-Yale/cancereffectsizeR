#' Set SBS signature weights
#' 
#' If you wish to use your own method to calculate sample-specific SBS signature weights
#' (as opposed the signature extraction built into trinuc_mutation_rates()), you can use
#' this function to load them into the CESAnalysis. Your input signatures will be used to
#' infer relative trinucleotide-context-specific mutation rates for all tumors. (This
#' means you can run set_signature_weights() or set_trinuc_rates(), but not both.) As in
#' trinuc_mutation_rates(), you can use a built-in set of signatures, such as COSMIC_v3.1,
#' or you can supply your own signature set definitions as documented in
#' ?trinuc_mutation_rates.
#'
#' The input data table must have a patient_id column and one column per
#' signature in the signature set. All samples in the CESAnalysis must be included in the
#' input table, and each sample's weights should have a sum on (0, 1]. Since these weights
#' are used by cancereffectsizeR to infer trinucleotide-context-specific relative rates
#' of SBS variants, each sample must have at least one non-artifact signature with
#' nonzero weight. (In the unlikely event that this is a problem, consider assigning
#' group-average signature weights to the artifact-only samples.)
#' 
#' @param cesa CESAnalysis
#' @param signature_set signature set name (see \code{list_ces_signature_sets()}), or a
#'   custom signature set (see documentation in \code{trinuc_mutation_rates()})
#' @param weights data.table of relative signature weights for each sample (see details)
#' @param ignore_extra_samples skip samples in the input table that are not in the
#'   CESAnalysis (when false, will stop with an error)
#' @export
set_signature_weights = function(cesa, signature_set, weights, ignore_extra_samples = FALSE) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be CESAnalysis.")
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
  
  previous_set = cesa@advanced$sbs_signatures[[signature_set_name]]
  if (! is.null(previous_set)) {
    if(! all.equal(previous_set, signature_set_data, check.attributes = F)) {
      stop("A signature set with the same name has already been used in the CESAnalysis, but it is not ",
           "identical to the input set.")
    }
  }
  
  if(! is(weights, "data.table")) {
    stop("weights should be a data.table (you may need to convert with as.data.table()).")
  }
  if (! is.logical(ignore_extra_samples) || length(ignore_extra_samples) != 1) {
    stop("ignore_extra samples should be T/F")
  }
  
  if (uniqueN(colnames(weights)) != length(weights)) {
    stop("Input weights table has repeated column names.")
  }
  
  if (all(c("total_sbs", "sig_extraction_sbs", "group_avg_blended") %in% colnames(weights))) {
    message("It looks like the input table came from a previous CES run. Meta columns like total_sbs will be dropped and regenerated.")
    weights = weights[, -c("total_sbs", "sig_extraction_sbs", "group_avg_blended")]
  }
  
  # validate weights
  signature_names = rownames(signatures)
  if (! identical(sort(c(signature_names, "patient_id")), sort(colnames(weights)))) {
    stop("Column names of weights must exactly match all signature names in the signature set (with none missing), ",
         "plus a patient_id column.")
  }
  
  if (! all(weights[, sapply(.SD, is.numeric), .SDcols = signature_names])) {
    stop("Not all signature weight columns in input are type numeric.")
  }
  if (! is(weights$patient_id, "character")) {
    stop("Input weights column patient_id should be type character.")
  }
  
  if(any(duplicated(weights$patient_id))) {
    stop("Some samples appear more than once in the weights table.")
  }
  
  
  if (length(setdiff(cesa@samples$patient_id, weights$patient_id)) > 0) {
    stop("Not all samples in the CESAnalysis appear in the input weights table.")
  }
  
  if (length(setdiff(weights$patient_id, cesa@samples$patient_id)) > 0 & ! ignore_extra_samples) {
    stop("There are samples in the input weight table that are not part of the CESAnalysis. (Re-run with ",
         "ignore_extra_samples = TRUE to ignore these.)")
  }
  
  # subset to CESAnalysis samples
  weights = weights[patient_id %in% cesa@samples$patient_id]
  
  if(anyNA(weights)) {
    stop("There are NA values in the input weights table.")
  }
  
  # Require sum of weights <= 1, but allow some floating point imprecision
  weights_okay = apply(weights[, -"patient_id"], 1, function(x) { total = sum(x); total > 0 && total <= 1 + 100 * .Machine$double.eps})
  if (! all(weights_okay)) {
    bad_samples = weights[which(! weights_okay), patient_id]
    stop("Some samples have all-zero weights or weights summing greater than 1:\n", paste(bad_samples, collapse = ", "), ".")
  }

  
  total_sbs_counts = cesa@maf[variant_type == "sbs"][, .(total_sbs = .N), by = "patient_id"]
  weights[total_sbs_counts, total_sbs := total_sbs, on = "patient_id"]
  weights[is.na(total_sbs), total_sbs := 0]
  weights[, c("sig_extraction_sbs", "group_avg_blended") := list(NA_integer_, NA)]
  setcolorder(weights, c("patient_id", "total_sbs", "sig_extraction_sbs", "group_avg_blended"))
  
  # Produce relative rates of biological process signatures
  # It's up to the user to figure out how to deal with a sample that has entirely artifact weighting
  artifact_signatures = signature_metadata[Likely_Artifact == TRUE, Signature]
  bio_weights = artifact_account(weights, artifact_signatures = artifact_signatures, 
                                          signature_names = signature_names, fail_if_zeroed = TRUE)
  
  trinuc_rates = calculate_trinuc_rates(as.matrix(bio_weights[, ..signature_names]), as.matrix(signatures), 
                                        bio_weights$patient_id)
  cesa = copy_cesa(cesa)
  cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix = trinuc_rates
  cesa@trinucleotide_mutation_weights$raw_signature_weights = weights
  cesa@trinucleotide_mutation_weights$signature_weights_table_with_artifacts = data.table(NULL)
  cesa@trinucleotide_mutation_weights$signature_weight_table = bio_weights
  cesa@advanced$sbs_signatures[[signature_set_name]] = copy(signature_set_data)
  cesa@samples[, sig_analysis_grp := 0L]
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}

#' Calculate relative rates of biological mutational processes
#' 
#' Sets artifact signature weights to zero and normalizes so that biologically-associated
#' weights sum to (1 - unattributed proportion) in each sample.
#' 
#' @param weights data.table of signature weights (can have extra columns)
#' @param artifact_signatures vector of artifact signature names (or NULL)
#' @param signature_names names of signatures in weights (i.e., all column names)
#' @param fail_if_zeroed T/F on whether to exit if a tumor would have all-zero weights.
#' @keywords internal
artifact_account = function(weights, signature_names, artifact_signatures = NULL, fail_if_zeroed = FALSE) {
  bio_weights = copy(weights)
  if(! is.null(artifact_signatures) && length(artifact_signatures) > 0) {
    bio_weights[, (artifact_signatures) := 0]
  }
  bio_weights[, adjust := sum(.SD), .SDcols = signature_names, by = "patient_id"]
  zeroed_out_index = bio_weights[adjust == 0, which = T]
  if (length(zeroed_out_index) > 0) {
    if(fail_if_zeroed) {
      stop("Some tumor(s) have all-zero weights across non-artifact signatures:\n",
           paste(bio_weights[zeroed_out_index, patient_id], collapse = ", "), ".")
    }
    bio_weights[zeroed_out_index, adjust := 1] # don't alter all-zero rows
  }
  bio_weights[, (signature_names) := lapply(.SD, `/`, adjust), .SDcols = signature_names]
  bio_weights[, adjust := NULL]
  return(bio_weights)
}

#' Calculate trinuc rates
#' 
#' Used internally to calculate trinuc rates from signature weights
#' 
#' If any relative rate is less than 1e-9, we add the lowest above-threshold rate to all
#' rates and renormalize rates so that they sum to 1.
#' 
#' @param weights matrix of signature weights
#' @param signatures matrix of signatures
#' @param tumor_names names of tumors corresponding to rows of weights
#' @return matrix of trinuc rates where each row corresponds to a tumor
#' @keywords internal
calculate_trinuc_rates = function(weights, signatures, tumor_names) {
  
  # get trinuc rates and normalize to relative rates
  trinuc_rates = weights %*% signatures
  trinuc_rates = t(apply(trinuc_rates, 1, function(x) {
    total_trinuc_count = sum(x)
    if (total_trinuc_count > 0) {
      return(x/total_trinuc_count)
    } else {
      return(x)
    }
  }))
  rownames(trinuc_rates) = tumor_names
  
  # if any sample has a rate under the floor of 1e-9 in any context (possible with the
  # right combination of signatures), add the lowest above-threshold rate to all and
  # renormalize
  add_pseudo = function(x) {
    if(any(x < 1e-9) && ! all(x < 1e-9)) {
      next_lowest = min(x[x >= 1e-9])
      x = x + next_lowest
      x = x / sum(x)
    }
    return(x)
  }
  
  # Have to keep transposing after apply
  trinuc_rates = t(apply(trinuc_rates, 1, add_pseudo))
  return(trinuc_rates)
}

#' Get MutationalPatterns contributions matrix
#' 
#' Reformat a signature weights table from mutational signature analysis into the
#' contributions matrix required for MutationalPatterns functions, including
#' visualizations.
#' @param signature_weight_table As created by trinuc_mutation_rates(); typically accessed
#'   via (CESAnalysis)$mutational_signatures.
#' @export
convert_signature_weights_for_mp = function(signature_weight_table) {
  if(! is(signature_weight_table, "data.table")) {
    stop("signature_weight_table should be data.table")
  }
  if (! "patient_id" %in% names(signature_weight_table)) {
    stop("Didn't find patient_id column. The input should be a table generated by trinuc_mutation_rates().")
  }
  cols_to_keep = setdiff(names(signature_weight_table), c("patient_id", "total_sbs",
                                                          "sig_extraction_sbs", "group_avg_blended"))
  rn = signature_weight_table$patient_id
  signature_weight_table = signature_weight_table[, ..cols_to_keep]
  contrib_mat = as.matrix(signature_weight_table)
  rownames(contrib_mat) = rn
  return(t(contrib_mat))
}
