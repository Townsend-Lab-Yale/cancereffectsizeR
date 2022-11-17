#'   Attribute cancer effects to mutational signatures
#'
#'   Within patients and across the cohort, calculates mutational source probabilities and the share
#'   of cancer effects attributable to each source, where sources are biologically-associated
#'   mutational signatures. See \href{https://doi.org/10.1093/molbev/msac084}{Attribution of Cancer Origins to Endogenous, Exogenous, and Preventable Mutational Processes}
#'   for background and applications.
#'
#' @param cesa CESAnalysis with cancer effects calculated for the variants of interest.
#' @param effects A table of cancer effect estimates for a set of variants, as produced with
#'   \code{ces_variant()}. Different sets of variants (or parameter choices in the
#'   \code{ces_variant} run) will affect output.
#' @param samples Samples for which to calculate mutational sources and effect shares; defaults to all
#'   samples. Only included samples will inform the cohort-level output.
#' @return A nested list containing mutational source probabilities for each variant (per patient
#'   and averaged over patients with the variant) and cancer effect shares attributable to each
#'   signature (per patient and across the cohort).
mutational_signature_effects <- function(cesa = cesa, effects = NULL, samples = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis object")
  }
  
  num_signature_sets = length(cesa$reference_data$snv_signatures)
  if(num_signature_sets == 0) {
    stop("Input CESAnalysis contains no associated signature definitions.")
  } 
  if(num_signature_sets > 1) {
    msg = paste0("Unusual situation: Input CESAnalysis has more than one set of associated SNV mutational signature definitions. ",
                 "It should be possible to re-do the analysis using just one signature set.")
    stop(pretty_message(msg, emit = F))
  }
  signature_defs = cesa$reference_data$snv_signatures[[1]]$signatures
  signature_names = rownames(signature_defs)
  
  # Get data.table with signature weights (one row per tumor)
  signature_weights = get_signature_weights(cesa)
  if (is.null(signature_weights)) {
    stop("Can't run because no signature weights are associated with this analysis.\n",
         "If you have them, you can add them to the CESAnalysis with set_signature_weights().")
  }
  if(nrow(signature_weights) != nrow(cesa$samples)) {
    stop("Not all samples in CESAnalysis have assigned signature weights, so can't run.")
  }
  
  # To reduce number of columns in output tables, drop signatures where weights are always 0.
  not_zero = signature_weights[, (sapply(.SD, function(x) any(x != 0))), .SDcols = signature_names]
  not_zero = names(which(not_zero == T))
  signature_weights = signature_weights[, c("Unique_Patient_Identifier", ..not_zero)]
  setkey(signature_weights, "Unique_Patient_Identifier")
  signature_names = signature_names[signature_names %in% not_zero]
  signature_defs = signature_defs[not_zero, ]
  
  if(! is.data.table(effects) || effects[, .N] == 0) {
    stop("Expected effects to be data.table.")
  } else {
    required_cols = c('variant_id', 'selection_intensity', 'variant_type')
    cols_to_use = names(effects)[names(effects) %in% required_cols]
    if(length(cols_to_use) < 3) {
      stop("Expected effects input to have variant_id, variant_type, and selection_intensity columns.")
    }
    
    if(! is.character(effects$variant_id) || ! is.character(effects$variant_type)) {
      stop('In effects input, expected variant_id and variant_type to be character.')
    }
    if(! is.numeric(effects$selection_intensity)) {
      stop('In effects input, expected selection_intensity to be numeric.')
    }
    
    if(length(cols_to_use) > 3) {
      stop("Found repeated column names in effects input (check variant_id, variant_type, and selection_intensity).")
    }
    if(uniqueN(effects$variant_id) != effects[, .N]) {
      stop("Some variant_id appear more than once in the input effects table.")
    }
    
    # Ensure that we're only handling AACs and SNVs
    other_variant_type_index = effects[! variant_type %in% c("snv", "aac"), which = T]
    if (length(other_variant_type_index) > 0) {
      effects = effects[! other_variant_type_index, ]
      warning("Some variants in effects input are neither coding mutations nor noncoding SNVs; these are being skipped.")
      
      # Make sure filtering didn't remove all variants
      if(effects[, .N] == 0) {
        stop("No variants remain.")
      }
    }
    
    missing = setdiff(effects$variant_id,
                      c(cesa@mutations$amino_acid_change$aac_id, cesa@mutations$snv$snv_id))
    if(length(missing) > 0) {
      msg = paste0("Some variants in effects input are not present in the CESAnalysis variant annotations. ",
                   "(This shouldn't happen. Make sure the effects table wasn't accidentally altered.)")
      stop(pretty_message(msg, emit = F))
    }
  }
  trinuc_rates = melt(cesa$trinuc_rates, id.vars = 'Unique_Patient_Identifier', variable.name = 'trinuc_context')

  
  # Collect SNV IDs: just the variant ID for SNVs; for AACs, pull information from mutation annotations
  effects[variant_type == "snv", snv_ids := list(list(variant_id)), by = "variant_id"]
  effects[variant_type == "aac", snv_ids := cesa@mutations$amino_acid_change[variant_id, constituent_snvs, on = 'aac_id']]
  
  variant_sources = effects[, .(snv = unlist(snv_ids), selection_intensity), by = "variant_id"]
  variant_sources[, trinuc_context := cesa@mutations$snv[snv, trinuc_mut, on = 'snv_id']]
  
  # Ensure that all SNVs have trinuc context annotation
  bad_snv_index = variant_sources[is.na(trinuc_context), which = T]
  if (length(bad_snv_index) > 0) {
    missing = variant_sources[bad_snv_index, variant]
    for (variant in missing) {
      warning(sprintf("Variant %s has incomplete SNV annotations, so it was skipped.", variant))
    }
    variant_sources = variant_sources[! variant %in% missing]
  }
  # Subset to desired samples
  if(is.null(samples)) {
    samples = cesa$samples$Unique_Patient_Identifier
  } else {
    samples = select_samples(cesa = cesa, samples = samples)$Unique_Patient_Identifier
  }
  
  # Subset to included samples and get a table of variants, SIs, UPIs
  maf = cesa$maf[samples, on = 'Unique_Patient_Identifier', nomatch = NULL]
  variant_sources = merge.data.table(variant_sources, maf[, .(snv = variant_id, Unique_Patient_Identifier)], 
                                          by = 'snv', all = FALSE)
  variant_sources[trinuc_rates, trinuc_rate := value, on = c('Unique_Patient_Identifier', 'trinuc_context')]
  
  
  # For each row (tumor-variant combination) in variant_snv_pairing, calculate
  # P(signature produced variant) = (sig_weight * sig_trinuc_rate) / tumor_trinuc_rate.
  sig_weights_by_snv = signature_weights[variant_sources$Unique_Patient_Identifier, .SD, 
                                     .SDcols = signature_names, on = 'Unique_Patient_Identifier']
  sig_trinuc_rates = t(signature_defs)[variant_sources$trinuc_context, signature_names]
  sig_prob_by_snv = (sig_weights_by_snv * sig_trinuc_rates) / variant_sources$trinuc_rate
  
  # Results reported for each UPI by variant_id, which includes a mix of SNV/AAC. For a given AAC
  # that has multiple associated SNV, a patient shouldn't have more than one SNV within the AAC
  # because it should typically have been called called as a multinucleotide variant in
  # preload_maf(). But if it does occur, the variant/UPI pairings will be repeated.
  # In the future, upstream improvements will prevent AACs from being called in such patients. 
  variant_sources = variant_sources[, .(variant_id, Unique_Patient_Identifier, selection_intensity)]
  variant_sources[, (signature_names) := sig_prob_by_snv]
  
  # For each variant, get the average source contributions across all patients with the variant.
  cohort_sources = variant_sources[, lapply(.SD, mean), 
                                          by = 'variant_id', .SDcols = signature_names]
  cohort_sources[variant_sources, selection_intensity := selection_intensity, on = 'variant_id']
  setcolorder(cohort_sources, c('variant_id', 'selection_intensity'))
  
  
  # Effect share for a given variant and signature is P(signature is source) * selection_intensity.
  # When reported for each patient and the cohort, we normalize such that sum across signatures is 1.
  effect_shares = variant_sources[, .(Unique_Patient_Identifier, .SD * variant_sources$selection_intensity), .SDcols = signature_names]
  
  # For each patient, sum the effect shares across all variants in the patient, and normalize.
  patient_effect_shares = effect_shares[, lapply(.SD, sum), .SDcols = signature_names, by = 'Unique_Patient_Identifier']
  patient_effect_shares[, (signature_names) := .SD/rowSums(.SD), .SDcols = signature_names, by = 'Unique_Patient_Identifier']
  
  # For all included samples, sum all effect shares and normalize.
  total_effect = sum(variant_sources$selection_intensity) # equivalent to sum(effect_shares[, .SD, .SDcols = signature_names])
  cohort_effect_shares = effect_shares[, lapply(.SD, function(x) sum(x)/total_effect), .SDcols = signature_names]
  cohort_effect_shares = setNames(as.numeric(cohort_effect_shares), names(cohort_effect_shares))
  effect_shares = list(patient = patient_effect_shares, cohort = cohort_effect_shares)
  
  return(list(variant_sources = list(patient = variant_sources, cohort = cohort_sources),
              effect_shares = list(patient = patient_effect_shares, cohort = cohort_effect_shares)))
}
