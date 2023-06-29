#' Assess variant prevalence and coverage
#' 
#' Determine variant prevalence (and how many samples have sequencing coverage) across your MAF data,
#' or within different groups of samples.
#'
#' @param cesa CESAnalysis
#' @param variant_ids variant names ("KRAS G12C") or full variant IDs. If left empty, uses
#'   non-overlapping variants as returned by `select_variants()` with \code{min_freq = 1}.
#' @param by Optionally, a vector of one or more sample table columns. Variant prevalence
#'   and coverage data will be broken down by the groups defined by unique combinations of
#'   values in these columns.
#' @export
variant_counts = function(cesa, variant_ids = character(), by = character()) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be CESAnalysis")
  }
  
  by_cols = by
  if(! is(by_cols, "character")) {
    stop("by should be character vector of CESAnalysis sample table column names.")
  }
  by_cols = unique(by_cols)
  missing_cols = setdiff(by_cols, names(cesa@samples))
  if(length(missing_cols) > 0) {
    msg = paste0("Column names given with by are not present in CESAnalysis sample table: ", 
                 paste0(missing_cols, collapse = ', '))
    stop(pretty_message(msg, emit = F))
  }
  
  if("Unique_Patient_Identifier" %in% by_cols) {
    stop("It doesn't really make sense to have Unique_Patient_Identifier in by.")
  }
  
  # Get sample table with just by_cols, UPI, covered_regions.
  # Handle case of covered_regions being in by_cols by renaming.
  samples = copy(cesa@samples)
  if("covered_regions" %in% by_cols) {
    by_cols[by_cols == "covered_regions"] = "cr_copy_for_by"
    samples[, cr_copy_for_by := covered_regions]
  }
  samples_with_by_cols = samples[, .SD, .SDcols = c(by_cols, "Unique_Patient_Identifier", "covered_regions")]
  
  if(! is(variant_ids, "character")) {
    stop("variant_ids should be character")
  }
  if(length(variant_ids) == 0) {
    variants = cesa$variants
    noncoding_snv_id = variants[variant_type == 'snv', variant_id]
    aac_ids = variants[variant_type == 'aac', variant_id]
    snv_from_aac = cesa@mutations$aac_snv_key[aac_ids, .(aac_id, snv_id), on = 'aac_id']
  } else {
    variants = sort_and_validate_variant_ids(cesa, variant_ids)
    noncoding_snv_id = variants[['snv_id']]
    aac_ids = variants[['aac_id']]
    unannotated = c(setdiff(aac_ids, cesa@mutations$amino_acid_change$aac_id),
                    setdiff(noncoding_snv_id, cesa@mutations$snv$snv_id))
    if (length(unannotated) > 0) {
      pretty_message(paste0(length(unannotated), " input variants are being left out of output ",
                            "because they're not annotated in the analysis. (Their MAF prevalence ",
                            "is zero, but samples covering hasn't been calculated.)"))
      aac_ids = setdiff(aac_ids, unannotated)
      noncoding_snv_id = setdiff(noncoding_snv_id, unannotated)
    }
    snv_from_aac = cesa@mutations$aac_snv_key[aac_ids, .(aac_id, snv_id), on = 'aac_id']
  }
  
  # call internal function
  return(.variant_counts(cesa, samples_with_by_cols, snv_from_aac, noncoding_snv_id, by_cols))
}

#' Internal variant prevalence and coverage calculation
#' 
#' Called by variant_counts() (and select_variants()) with validated inputs.
#' 
#' @param cesa CESAnalysis
#' @param samples validated samples table
#' @param snv_from_aac data.table with columns aac_id, snv_id (validated and with annotations in CESAnalysis)
#' @param noncoding_snv_id vector of snv_ids to treat as noncoding variants
#' @param by_cols validated column names from sample table that are suitable to use for counting by.
#' @keywords internal
.variant_counts = function(cesa, samples, snv_from_aac, noncoding_snv_id, 
                           by_cols = character()) {
  
  maf = copy(cesa@maf)
  if(maf[, .N] == 0) {
    # Empty maf lacks column names (will probably change soon)
    maf = data.table(Unique_Patient_Identifier = character(), variant_id = character())
  }
  
  
  # Get counts by by_cols, with 0-count combinations also listed
  get_complete_counts = function(dt) {
    final_counts = setDT(samples[, expand.grid(c(list(variant_id = unique(dt$variant_id)), 
                                             lapply(.SD, unique)), stringsAsFactors = FALSE), 
                             .SDcols = by_cols])
    final_counts[, N := 0]
    dt = merge.data.table(dt, samples[, .SD, .SDcols = c(by_cols, 'Unique_Patient_Identifier')], all.x = TRUE, by = 'Unique_Patient_Identifier')
    nonzero_counts = dt[! is.na(Unique_Patient_Identifier), .(count = .N), by = c('variant_id', by_cols)]
    final_counts[nonzero_counts, N := count, on = c(by_cols, 'variant_id')]
    return(final_counts)
  }
  
  
  combined_counts = data.table()
  if(snv_from_aac[, .N] > 0) {
    snv_from_aac_counts = setDT(maf[snv_from_aac, .(Unique_Patient_Identifier, variant_id = aac_id), 
                                    on = c(variant_id = "snv_id"), 
                                    allow.cartesian = T])
    aac_count_output = get_complete_counts(dt = snv_from_aac_counts)
    combined_counts = rbind(combined_counts, aac_count_output)
  }
  
  
  # Get prevalences of SNVs
  if(length(noncoding_snv_id) > 0) {
    snv_counts = maf[noncoding_snv_id, .(Unique_Patient_Identifier, variant_id), on = 'variant_id']
    final_snv_counts = get_complete_counts(snv_counts)
    combined_counts = rbind(combined_counts, final_snv_counts)
  }
  
  
  if(combined_counts[, .N] == 0) {
    stop("No variants to count.")
  }
  
  calc_cov = function(variant_id) {
    # Count samples with coverage in each group
    covering_by_variant = mget(variant_id, cesa@mutations$variants_to_cov)
    
    dts = list(copy(samples))
    if(length(by_cols) > 0) {
      dts = split(samples, by = by_cols)
    }
    
    cov_by_group = rbindlist(lapply(dts, function(dt) {
      # Hash sample counts for each panel for the specific by-group
      sample_count_by_cr = dt[, .N, by = 'covered_regions']
      num_generic_wg = sample_count_by_cr[covered_regions == 'genome', N] # to add in at the end

      if(length(num_generic_wg) == 0) {
        num_generic_wg = 0
      }
      sample_count_by_cr = list2env(as.list(setNames(sample_count_by_cr$N, 
                                                     nm = sample_count_by_cr$covered_regions)),
                                    parent = emptyenv())
      cov_by_group = lapply(covering_by_variant, 
                            function(x) sum(unlist(mget(x, sample_count_by_cr, ifnotfound = 0))))
      cov_by_group = data.table(variant_id = names(cov_by_group), num_cov = unlist(cov_by_group))
      sapply(by_cols, function(x) cov_by_group[, (x) := dt[[x]][1]])
      
      # Add in the generic whole-genome samples
      cov_by_group[, num_cov := num_cov + num_generic_wg]
      return(cov_by_group)
    }))
    return(cov_by_group)
  }
  combined_cov_output = calc_cov(unique(combined_counts$variant_id))
  
  full_output = merge.data.table(combined_counts, combined_cov_output, by = c('variant_id', by_cols), sort = F)
  return(full_output[order(variant_id)])
}

#' Sort and validate input variant IDs
#'
#' Sorts input variant IDs by type, completes IDs by adding protein ID to plain variant
#' names (e.g. "KRAS G12C"), and ensures that IDs are valid even if not present in
#' annotations. This includes verifying that reference alleles are correct in SNV IDs and
#' that amino-acid-changes are possible given the coding sequence.
#' 
#' @return List of variant_ids, with each element corresponding to one variant_type.
#' @param cesa CESAnalysis
#' @param drop_unannotated Whether to include variants that are not annotated in output.
#' @param input_ids Variant names/IDs, typically from user.
#' @keywords internal
sort_and_validate_variant_ids = function(cesa, input_ids, drop_unannotated = FALSE) {
  if(! is.character(input_ids)) {
    stop("input_ids should be character vector of variant identifiers")
  }
  input_ids = unique(na.omit(input_ids))
  if(length(input_ids) == 0) {
    stop("input_ids doesn't list any variants.")
  }
  snv_ids = intersect(input_ids, cesa@mutations$snv$snv_id)
  input_ids = setdiff(input_ids, snv_ids)
  
  aac_ids = intersect(input_ids, cesa@mutations$amino_acid_change$aac_id)
  input_ids = setdiff(input_ids, aac_ids)
  
  ## To-do: insert DBS/indel logic here (and below)
  dbs_ids = intersect(input_ids, cesa@mutations$dbs$dbs_id)
  input_ids = setdiff(input_ids, dbs_ids)
  
  dbs_aac_ids = intersect(input_ids, cesa@mutations$dbs_codon_change$dbs_aac_id)
  input_ids = setdiff(input_ids, dbs_aac_ids)
  
  
  if(length(input_ids) > 0) {
    input_ids = sub(' ', '_', input_ids)
    tmp_dt = copy(cesa@mutations$amino_acid_change)
    tmp_dt[, tmp_name :=  paste(gene, aachange, sep = '_')]
    more_aac_ids = tmp_dt[input_ids, aac_id, on = 'tmp_name', nomatch = NULL]
    input_ids = setdiff(input_ids, tmp_dt[more_aac_ids, tmp_name, on = 'aac_id'])
    aac_ids = union(aac_ids, more_aac_ids)
  }
  
  # Keep going if all unmatched IDs look valid. (Presumably, these are absent from annotations.)
  if(length(input_ids) > 0) {
    apparent_snvs = input_ids[grepl(':\\d+_[ACTG]>[ACGT]$', input_ids)]
    if(length(apparent_snvs) > 0) {
      validate_snv_ids(apparent_snvs, get_cesa_bsg(cesa))
    }
    if(! drop_unannotated) {
      snv_ids = union(snv_ids, apparent_snvs)
    }
    # insert more indel/dbs logic
    apparent_aac = setdiff(input_ids, apparent_snvs)
    apparent_aac = complete_aac_ids(apparent_aac, .ces_ref_data[[cesa@ref_key]])
    aac_problems = validate_aac_ids(apparent_aac, .ces_ref_data[[cesa@ref_key]])
    if(! is.null(aac_problems)) {
      print(aac_problems)
      stop("The above variant IDs from input_ids appear invalid.")
    }
    if(! drop_unannotated) {
      aac_ids = union(aac_ids, apparent_aac)
    }
  }
  
  return(list(snv_id = snv_ids, aac_id = aac_ids))
}

