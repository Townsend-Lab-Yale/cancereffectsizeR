#' Create CompoundVariantSet directly from SNV IDs
#'
#' A compound variant is an arbitrary group of SNVs that can be tested for selection as if
#' it were a single variant. (Any sample with one or more of the constituent SNVs "has the
#' compound variant", and the baseline rate of the variant is the rate of having at least
#' one of the SNVs.) A CompoundVariantSet is a collection of disjoint compound variants;
#' that is, no SNV can appear in more than one compound variant.\cr
#' Example: \code{CompoundVariantSet(cesa, variant_id = list(kras12 = c("KRAS G12C", "KRAS G12D",
#' "KRAS G12V")))} creates a CompoundVariantSet containing one compound variant. To create
#' a large set, it's usually easier to use define_compound_variants(), which calls this
#' function internally to define compound variants from an input variant table.
#' 
#' If you're using this function because you have a complex use case that
#' define_compound_variants() can't handle, please let us know so we can try to make
#' improvements!
#' 
#' @param cesa CESAnalysis (used to access variant annotations)
#' @param variant_id Vector of variant IDs to include in one compound variant, or a list of
#'   vectors, each of which defines a separate compound variant. If the vector or list is
#'   named, names will be kept. Otherwise, compound variants will be named sequentially.
#' @export
CompoundVariantSet = function(cesa, variant_id) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  
  if (cesa@mutations$snv[, .N] == 0) {
    stop("There are no mutation annotations in the CESAnalysis")
  }
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  if (is(variant_id, "character")) {
    prev_name = names(variant_id)
    variant_id = as.list(unname(variant_id))
    names(variant_id) = prev_name
  }
  if(! is(variant_id, "list")) {
    stop("variant_id should be list or character")
  }
  
  if(length(variant_id) == 0) {
    stop("variant_id can't be 0-length.")
  }
  if (! all(sapply(variant_id, is.character))) {
    stop("All elements of variant_id must be type character")
  }
  if (! all(sapply(variant_id, length) > 0)) {
    stop("All elements of snv_id list must have nonzero length")
  }
  
  compound_names = names(variant_id)
  legal_compound_name = '^[a-z0-9][-0-9a-z\\_\\.]*$'
  if (! is.null(compound_names) && ! all(grepl(legal_compound_name, compound_names, ignore.case = T))) {
    stop("Some or all compound variant names are invalid. Must start with letter or number and contain ",
         "only letters, numbers, hyphens, underscores, periods.")
  }
  
  # first check that no variants overlap (will need to check again after breaking down AACs)
  all_ids = unlist(variant_id)
  if (length(unique(all_ids)) != length(all_ids)) {
    dup = which(duplicated(all_ids))[1]
    stop("Some SNVs appear multiple times in the input set (e.g., ", dup, "). If you want to use overlapping compound variants, ",
         "call create_compound_variants() multiple times to make separate variant sets.")
  }
  
  
  # replace spaces with underscores (mainly to allow things like KRAS G12C -> KRAS_G12C)
  variant_id = lapply(variant_id, function(x) gsub(" ", "_", x))
  
  selected = select_variants(cesa, variant_ids = all_ids, include_subvariants = T)
  selected_snvs = selected[variant_type == "snv"]
  setkey(selected_snvs, "variant_id")
  selected_aacs = selected[variant_type == "aac"]
  setkey(selected_aacs, "variant_id")
  setindex(selected_aacs, "variant_name")
  
  # count how many compound variants have no covered_regions coverage,
  # and also have how many have coverage in only some covered_regions
  no_shared_coverage_count = 0
  some_incomplete_coverage = 0
  
  # count how many distinct covered_region sets have samples (excluding the whole-genome-coverage set)
  num_coverages_with_sample = cesa@samples[covered_regions != "genome", length(unique(covered_regions))]
  compound_variants = list()
  for (i in 1:length(variant_id)) {
    current_ids = variant_id[[i]]
    current_snvs = selected_snvs[current_ids, nomatch = NULL]
    current_snv_ids = current_snvs$variant_id
    current_coverage = mget(current_snv_ids, envir = cesa@mutations$variants_to_cov)
    if (current_snvs[, .N] < length(current_ids)) {
      remaining_ids = setdiff(current_ids, current_snvs$variant_id)
      current_aacs = selected_aacs[remaining_ids, nomatch = NULL]
      remaining_ids = setdiff(remaining_ids, current_aacs$variant_id)
      
      # handle user input of variant names (in place of IDs)
      if (length(remaining_ids) > 0) {
        current_aacs = rbind(current_aacs, selected_aacs[remaining_ids, nomatch = NULL, on = "variant_name"])
      }
      
      if (current_aacs[, .N] + current_snvs[, .N] != length(current_ids)) {
        stop("Problem finding unique and complete annotations for variant ", i, " in set.")
      }
      snv_ids_from_aac = current_aacs[, unlist(constituent_snvs)]
      if (any(duplicated(snv_ids_from_aac))) {
        stop("Item ", i, " of input contains AAC variants with overlapping constituent SNVs.")
      }
      current_snv_ids = c(current_snv_ids, snv_ids_from_aac)
      if (any(duplicated(current_snv_ids))) {
        stop("Item ", i, " of input has overlapping SNV/AAC IDs.")
      }
      current_coverage = c(current_coverage, mget(current_aacs$variant_id, 
                                                  envir = cesa@mutations$variants_to_cov))
    }
    minimal_coverage = Reduce(intersect, current_coverage)
    maximal_coverage_size = length(Reduce(union, current_coverage))
    min_cov_size = length(minimal_coverage)
    if (min_cov_size == 0) {
      no_shared_coverage_count = no_shared_coverage_count + 1
    }
    if (min_cov_size < maximal_coverage_size) {
      some_incomplete_coverage = some_incomplete_coverage + 1
    }
    compound_variants[[i]] = list(snv_id = current_snv_ids, coverage = minimal_coverage)
    variant_id[[i]] = current_snv_ids
  }
  
  # first check that no variants overlap (will need to check again after breaking down AACs)
  if (any(duplicated(unlist(variant_id)))) {
    stop("After breaking down AACs into SNVs, some compound variants share overlapping SNVs. ",
         "If this is desired, use create_compound_variants() multiple times to make separate variant sets.")
  }
  
  # name sequentially unless user provided names
  if (is.null(names(variant_id))) {
    names(compound_variants) = paste0("comp.", 1:length(variant_id))
  } else {
    names(compound_variants) = names(variant_id)
  }
  
  set_size = length(variant_id)
  if (some_incomplete_coverage > 0) {
    msg = paste0("\nNote: ", some_incomplete_coverage, " of ", set_size, " compound variants ",
                 "have sequencing coverage in different sets of samples at different constituent variant sites. ",
                 "Selection functions will only use samples with coverage at all sites. If you want to ",
                 "include more samples, remove uncovered variants from these compound variants.")
    pretty_message(msg, black = F)
  }
  
  if(no_shared_coverage_count > 0 && cesa@samples[covered_regions == "genome", .N] == 0) {
    warning(no_shared_coverage_count, " of ", set_size, " compound variants are not fully covered ",
            "by any sample's covered_regions. Remove uncovered variants to be able to test these for selection.")
  }
  
  # convert list to a data.table for easier handling
  variants = lapply(compound_variants, "[[", 1)
  coverage = lapply(compound_variants, "[[", 2)
  
  compound_snvs = data.table(compound_name = names(compound_variants), variants = variants)
  compound_snvs = compound_snvs[, .(snv_id = unlist(variants)), by = "compound_name"]
  compound_snvs[, shared_cov := coverage[compound_name]]
  
  # add in simplified AAC/gene annotations
  aac_anno = cesa@mutations$aac_snv_key[compound_snvs$snv_id, on = 'snv_id', nomatch = NULL]
  aac_anno[cesa@mutations$amino_acid_change, c("gene", "aachange") := list(gene, aachange), on = 'aac_id']
  aac_anno[, variant_name := paste0(gene, '_', aachange)]
  aac_anno = aac_anno[, .(variant_names = list(unique(variant_name)), genes = list(unique(gene)), 
                          aachanges = list(unique(aachange)), num_genes = uniqueN(gene)), by = 'snv_id']
  aac_anno[num_genes > 1, aachanges := variant_names]
  compound_snvs[aac_anno, c('genes', 'aachanges') := list(genes, aachanges), on = 'snv_id']

  sample_table = cesa@maf[compound_snvs$snv_id, Unique_Patient_Identifier, on = "variant_id", by = "variant_id", nomatch = NULL]
  sample_table[cesa@samples, covered_regions := covered_regions, on = "Unique_Patient_Identifier"]
  sample_table[compound_snvs, c("compound_name", "shared_cov") := list(compound_name, shared_cov), on = c(variant_id = "snv_id")]
  sample_table[, compound_covered := covered_regions == "genome" | covered_regions %in% unlist(shared_cov), by = "variant_id"]
  total_freq = sample_table[, .(total_maf_freq = .N), by = "variant_id"]
  
  covered_sample_table = sample_table[compound_covered == T]
  covered_freq = covered_sample_table[, .(shared_cov_maf_freq = .N), by = "variant_id"]
  
  # handle frequencies
  compound_snvs[covered_freq, shared_cov_maf_freq := shared_cov_maf_freq, on = c(snv_id = "variant_id")]
  compound_snvs[is.na(shared_cov_maf_freq), shared_cov_maf_freq := 0]
  compound_snvs[total_freq, total_maf_freq := total_maf_freq, on = c(snv_id = "variant_id")]
  compound_snvs[is.na(total_maf_freq), total_maf_freq := 0]
  
  
  # also create a table summarizing each compound variant; also drop shared_cov from compound_snvs
  compound_stats = compound_snvs[, .(num_snv = .N, shared_cov = shared_cov[1], shared_cov_subvariant_freq = sum(shared_cov_maf_freq),
                                    total_subvariant_freq = sum(total_maf_freq)), by = "compound_name"]
  compound_snvs[, shared_cov := NULL]
  
  # shared_cov_freq is the number of samples that have one or more of a compound variant, within samples have coverage
  # at all compound variant sites
  compound_counts_cov = unique(covered_sample_table, by = c("variant_id", "Unique_Patient_Identifier"))
  compound_counts_cov = compound_counts_cov[, .(shared_cov_freq = uniqueN(Unique_Patient_Identifier)), by = "compound_name"]
  compound_stats = compound_stats[compound_counts_cov, on = "compound_name"]
  compound_counts_total = unique(sample_table, by = c("variant_id", "Unique_Patient_Identifier"))
  
  # Record which samples (including those outside shared coverage) have the variant
  samples_by_variant = compound_counts_total[, .(samples = list(unique(Unique_Patient_Identifier))), by = "compound_name"]
  samples_with = stats::setNames(samples_by_variant$samples, samples_by_variant$compound_name)
  compound_counts_total = samples_by_variant[, .(total_freq = length(samples[[1]])), by = "compound_name"]
  compound_stats = compound_stats[compound_counts_total, on = "compound_name"]
  setcolorder(compound_stats, c("compound_name", "num_snv", "shared_cov", "shared_cov_freq", "total_freq"))

  

  # For simplicity, we don't want this variant set to be used after adding new samples,
  # even if there aren't new covered_regions
  num_samples = cesa@samples[, .N]
  
  setkey(compound_stats, "compound_name")
  setkey(compound_snvs, "compound_name", physical = FALSE)
  comp = new("CompoundVariantSet", snvs = compound_snvs, compounds = compound_stats, sample_calls = samples_with,
             cesa_uid = cesa@advanced$uid, cesa_num_samples = num_samples)
  return(comp[compound_names]) # Make order match input
}


#' Divide batches of variants into a CompoundVariantSet
#'
#' A compound variant is an arbitrary group of SNVs that can be tested for selection as if
#' it were a single variant. (Any sample with one or more of the constituent SNVs "has the
#' compound variant", and the baseline rate of the variant is the rate of having at least
#' one of the SNVs.) A CompoundVariantSet is a collection of disjoint compound variants;
#' that is, no SNV can appear in more than one compound variant. After collecting variants
#' of interest into a table using select_variants()--and further subsetting or annotating
#' the table as desired--use this function to produce a CompoundVariantSet that combines
#' variants into distinct compound variants based on your criteria.
#' 
#' This function works first by splitting the input table by the columns given in
#' \code{by}. For example, splitting on "gene" will split the table into gene-specific
#' subtables. Then, each subtable is divided into compound variants based on
#' \code{merge_distance}. All variants in each subtable within the specified genomic
#' distance of each other will be merged into a candidate compound variant, and then
#' compound variants will be repeatedly merged until the nearest two variants in each pair
#' of compound variants are not within \code{merge_distance}. Note that overlapping
#' variants will always be merged unless you use \code{by} to separate them into different
#' subtables (for example, by splitting on alt or aa_alt). If you use \code{by} to split
#' variants by some functional annotation, you can set \code{merge_distance} very high to
#' merge all same-chromosome sites (e.g., 1e9 on human genome). To merge sites across chromosomes,
#' set \code{merge_distance = Inf}.
#' 
#' @param cesa CESAnalysis
#' @param variant_table Data table of variants, in the style generated by select_variants().
#' @param by One or more column names to use for initial splitting of the input table into variant
#' groups. Each distinct group will then be further divided into compound variants based on \code{merge_distance}
#' @param merge_distance maximum genomic distance between a given variant and the nearest
#'   variant in compound variant for the variant to variant to be merged into the compound
#'   variant (as opposed to being assigned to its own compound variant).
#' @export
define_compound_variants = function(cesa, variant_table, by = NULL, merge_distance = 0) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis.")
  }
  if (! is(variant_table, "data.table")) {
    stop("variant_table should be a data.table.")
  }
  
  if (! is.numeric(merge_distance) || length(merge_distance) != 1) {
    stop("merge_distance should be positive integer (or Inf).")
  }
  if(! merge_distance == Inf) {
    if(as.integer(merge_distance) - merge_distance != 0 || merge_distance < 0) {
      stop("merge_distance should be positive integer (or Inf).")
    }
    merge_distance = as.integer(merge_distance)
  }
  
  
  if (! all(c("chr", "start", "end", "center_nt_pos") %in% names(variant_table))) {
    stop("variant_table should have chr/start/end/center_nt_pos (as from select_variants()).")
  }
  
  by_cols = by
  if(! is.null(by_cols)) {
    if (! is(by_cols, "character")) {
      stop("by should be character (must be column names of variant_table).")
    }
    if (! all(by_cols %in% names(variant_table))) {
      stop("not all column names specified with \"by\" are present in variant_table.")
    }
    
    # deal with NAs in "by" columns by converting to character and giving informative name
    if (variant_table[, anyNA(.SD), .SDcols = by_cols]) {
      variant_table = copy(variant_table)
      for (col in by_cols) {
        variant_table[, (col) := as.character(variant_table[[col]])]
        variant_table[is.na(variant_table[[col]]), (col) := paste0(col, ".NA")]
      }
    }
    
    table_list = split(variant_table, by = by_cols)
  } else {
    table_list = list(comp = variant_table)
  }
  
  variant_chunks = list()
  num_split_groups = length(table_list)
  split_group_names = names(table_list)
  
  for (i in 1:num_split_groups) {
    variant_table = table_list[[i]]
    curr_split_group_name = split_group_names[[i]]
    # If merge_distance == Inf, combine all variants in group
    if (merge_distance == Inf || variant_table[, .N] == 1) {
      current_group = list(variant_table$variant_id)
      names(current_group) = paste0(curr_split_group_name, '.1')
      variant_chunks = c(variant_chunks, current_group)
      next
    }
    non_aac = variant_table[variant_type != 'aac', .(variant_id, chr, start, end)]
    aac_table = variant_table[variant_type == 'aac']
    aac_table = rbindlist(list(aac_table[, .(variant_id, chr, start, end = start)], 
                               aac_table[, .(variant_id, chr, start = end, end)],
                               aac_table[, .(variant_id, chr, start = center_nt_pos, end = center_nt_pos)]))
    sites_by_variant = rbind(non_aac, aac_table)
    
    
    # Use reduce to merge variants that overlap the same sites (within user-supplied merge_distance)
    # Since some variants may only overlap partially, we'll get a list of character vectors
    # in which some variants appear in more than one vector. From there, we go through
    # all elements of the list merging each pair of subgroups that has a shared variant,
    # continuing to cycle through until no more groups need to be merged.
    gr = GenomicRanges::makeGRangesFromDataFrame(sites_by_variant, keep.extra.columns = T)
    reduced = GenomicRanges::reduce(gr, with.revmap = T, min.gapwidth = merge_distance)
    curr_groups = unique(lapply(reduced$revmap, function(x) unique(sites_by_variant$variant_id[x])))
    
    batch_num = new.env() # hash variant ID to "batch number" (which element of list variant list it belongs in)
    next_batch_num = 1
    while (TRUE) {
      any_changed = FALSE
      for (group in curr_groups) {
        for (variant in group) {
          # check if each variant already has a place assigned in list
          previous_batch_num = batch_num[[variant]]
          if (! is.null(previous_batch_num)) {
            break
          }
        }
        if(is.null(previous_batch_num)) {
          for (variant in group) {
            batch_num[[variant]] = next_batch_num
          }
          next_batch_num = next_batch_num + 1
        } else {
          for (variant in group) {
            curr_var_batch = batch_num[[variant]]
            if(is.null(curr_var_batch) || curr_var_batch != previous_batch_num) {
              batch_num[[variant]] = previous_batch_num
              any_changed = T
            }
          }
        }
      }
      if (! any_changed) break;
    }
    merged_groups = vector(mode = "list", length = next_batch_num - 1)
    for (variant in ls(batch_num)) {
      num = batch_num[[variant]]
      merged_groups[[num]] = c(merged_groups[[num]], variant)
    }
    names(merged_groups) = paste(curr_split_group_name, 1:length(merged_groups), sep = ".")
    variant_chunks = c(variant_chunks, merged_groups)
  }
  return(CompoundVariantSet(cesa, variant_id = variant_chunks))
}