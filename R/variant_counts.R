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
  
  # get sample table with just by_cols, UPI, covered_regions
  # handle case of covered_regions being in by_cols by renaming
  # And yes, convert by_cols to factor to more easily complete output table
  samples = copy(cesa@samples)
  if("covered_regions" %in% by_cols) {
    by_cols[by_cols == "covered_regions"] = "cr_copy_for_by"
    samples[, cr_copy_for_by := covered_regions]
  }
  samples_with_by_cols = samples[, .SD, .SDcols = c(by_cols, "Unique_Patient_Identifier", "covered_regions")]
  if(length(by_cols) > 0) {
    samples_with_by_cols[, (by_cols) := lapply(.SD, as.factor), .SDcols = by_cols]
  }
  
  # deal with data.table dcast bug (besides, "NA" and <NA> would conflict anyway in output headers)
  for (col in by_cols) {
    if (anyNA(samples_with_by_cols[[col]])) {
      if ("NA" %in% levels(samples_with_by_cols[[col]])) {
        stop("In column ", col, ", there are both missing values (NAs) and entries matching \"NA\".")
      }
      else {
        levels(samples_with_by_cols[[col]]) = c(levels(samples_with_by_cols[[col]]), "NA")
        samples_with_by_cols[is.na(samples_with_by_cols[[col]]), (col) := "NA"]
      }
    }
  }
  
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
  
  # Get prevalences of AACs
  aac_counts = data.table()
  aac_count_output = data.table()
  
  if(snv_from_aac[, .N] > 0) {
    snv_from_aac_counts = setDT(maf[snv_from_aac, .(Unique_Patient_Identifier, aac_id), on = c(variant_id = "snv_id"), nomatch = NULL])
    snv_from_aac_counts = merge.data.table(snv_from_aac_counts, samples, by = "Unique_Patient_Identifier", all.x = TRUE)
    aac_counts = snv_from_aac_counts[, .N, by = c(by_cols, 'aac_id')]
    setnames(aac_counts, 'aac_id', 'variant_id')
    zero_prev_aac = setdiff(snv_from_aac$aac_id, aac_counts$variant_id)
    if (length(zero_prev_aac) > 0) {
      aac_counts = rbind(aac_counts, data.table(variant_id = zero_prev_aac, N = 0), fill = T)
    }
    if(length(by_cols) > 0) {
      aac_count_output = dcast.data.table(aac_counts, variant_id ~ ..., value.var = "N", fill = 0, drop = F)
    } else {
      aac_count_output = copy(aac_counts)
    }
  }
  
  
  # Get prevalences of SNVs
  combined_count_output = aac_count_output
  snv_counts = data.table()
  snv_count_output = data.table()
  if(length(noncoding_snv_id) > 0) {
    snv_counts = maf[noncoding_snv_id, .(Unique_Patient_Identifier, variant_id), on = 'variant_id', nomatch = NULL]
    snv_counts = merge.data.table(snv_counts, samples, by = "Unique_Patient_Identifier", all.x = TRUE, sort = F)
    snv_counts = snv_counts[, .N, by = c(by_cols, 'variant_id')]
    zero_prev_snvs = setdiff(noncoding_snv_id, snv_counts$variant_id)
    if (length(zero_prev_snvs) > 0) {
      snv_counts = rbind(snv_counts, data.table(variant_id = zero_prev_snvs, N = 0), fill = T)
    }
    if(length(by_cols) == 0) {
      snv_count_output = copy(snv_counts)
    } else {
      snv_count_output = dcast.data.table(snv_counts, variant_id ~ ..., value.var = "N", fill = 0, drop = F)
    }
    combined_count_output = rbind(aac_count_output, snv_count_output)
  }
  
  if(combined_count_output[, .N] == 0) {
    stop("No variants to count.")
  }
  prevalence_cols = setdiff(names(combined_count_output), 'variant_id')
  new_prevalence_cols = 'total_prevalence'
  if(length(by_cols) > 0) {
    new_prevalence_cols = paste0(prevalence_cols, '_prevalence')
  }
  setnames(combined_count_output, prevalence_cols, new_prevalence_cols)
  
  
  calc_cov = function(count_output, annotations) {
    # Count samples with coverage in each group
    # Start by cross-joining all by_col factor combinations with all variant IDs (as given in count output)
    full_cov = do.call(CJ, c(as.list(samples[, .SD, .SDcols = by_cols]), 
                             list(variant_id = count_output$variant_id), unique = T))
    
    # stringi::stri_join_list can't handle character(0), so change to '',
    # include in 1-length annotations edge case
    if(annotations[, .N] == 1) {
      if(length(annotations$covered_in[[1]]) == 0) {
        annotations$covered_in[[1]] = list('')
      }
    } else {
      annotations[sapply(covered_in, length) == 0, covered_in := ''] 
    }
    
    annotations[, flat_cov := stringi::stri_join_list(covered_in)]
    full_cov[annotations, flat_cov := flat_cov, on = 'variant_id']
    
    # We only need to find set of coverage counts once per unique covered_in combination,
    # and then can copy the counts into full_cov
    small_cov = full_cov[ , .SD[1], by = c(by_cols, "flat_cov")]
    small_cov[annotations, covered_in := covered_in, on = 'variant_id'] # much smaller join
    cov_counts_by_group = samples[, .N, by = c(by_cols, "covered_regions")]
    setkeyv(cov_counts_by_group, c(by_cols, "covered_regions"))
    small_cov[, rn := 1:.N]
    
    # break open covered_in lists, making one row per variant/covering-region pair
    by_cov_region = small_cov[, .(covered_regions = c(unlist(covered_in), "genome")), by = 'rn']
    by_cov_region = merge.data.table(by_cov_region, small_cov[, -"covered_in"], by = 'rn', sort = F)
    setkeyv(by_cov_region, c(by_cols, "covered_regions"))
    by_cov_region[, num_cov := cov_counts_by_group[by_cov_region, N]]
    cov_counts = by_cov_region[, .(num_cov = sum(num_cov, na.rm = T)), by = "rn"]
    small_cov[cov_counts, num_cov := num_cov, on = 'rn']
    full_cov[small_cov, num_cov := num_cov, on = c(by_cols, 'flat_cov')]
    
    if(length(by_cols) > 0) {
      by_site = full_cov[, .(total_covering = sum(num_cov)), by = 'variant_id']
      cov_output = dcast.data.table(full_cov[, -c("flat_cov")], 
                                    variant_id ~ ..., value.var = "num_cov", fill = 0, drop = F)
      cov_output[by_site, total_covering := total_covering, on = 'variant_id']
    } else {
      cov_output = full_cov[, .(variant_id, num_cov)]
    }
    return(cov_output)
  }
  
  # Copy needed part of annotations
  
  aac_anno = cesa@mutations$amino_acid_change[snv_from_aac[, unique(aac_id)]]
  snv_anno = cesa@mutations$snv[noncoding_snv_id, ]
  setnames(aac_anno, 'aac_id', 'variant_id')
  setnames(snv_anno, 'snv_id', 'variant_id')
  
  combined_cov_output = data.table()
  if(aac_count_output[, .N] > 0) {
    aac_cov_output = calc_cov(aac_count_output, aac_anno)
    aac_cov_output[, variant_type := 'aac']
    combined_cov_output = aac_cov_output
  }
  if(snv_count_output[, .N] > 0) {
    snv_cov_output = calc_cov(snv_count_output, snv_anno)
    snv_cov_output[, variant_type := 'snv']
    combined_cov_output = rbind(combined_cov_output, snv_cov_output)
  }
  
  cov_cols = setdiff(names(combined_cov_output), c('variant_id', 'variant_type', 'total_covering'))
  new_cov_cols = paste0(cov_cols, '_covering')
  if(length(by_cols) == 0) {
    new_cov_cols = 'total_covering'
  }
  setnames(combined_cov_output, cov_cols, new_cov_cols)
  
  # Add total counts for both, if using by_col (already happened for coverage in helper function)
  first_cols = c('variant_id', 'variant_type')
  if (length(by_cols) > 0) {
    site_snv_counts = data.table()
    if(snv_counts[, .N] > 0) {
      site_snv_counts = snv_counts[, .(total_prevalence = sum(N)), by = 'variant_id']
    }
    site_aac_counts = data.table()
    if(aac_counts[, .N] > 0) {
      site_aac_counts = aac_counts[, .(total_prevalence = sum(N)), by = 'variant_id'] 
    }
    combined_count_output[rbind(site_snv_counts, site_aac_counts), total_prevalence := total_prevalence, on = 'variant_id']
    first_cols = c(first_cols, 'total_prevalence', 'total_covering')
  }
  
  full_output = merge.data.table(combined_count_output, combined_cov_output, by = 'variant_id', sort = F)
  setcolorder(full_output, c(first_cols, unlist(S4Vectors::zipup(new_prevalence_cols, new_cov_cols))))
  return(full_output)
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

