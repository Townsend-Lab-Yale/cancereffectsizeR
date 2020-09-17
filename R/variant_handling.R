#' Select and filter variants
#'
#' This is a lightweight function to help you find and manipulate variant data from your
#' CESAnalysis's MAF data and mutation annotation tables. This function has a variant
#' gathering step followed by a filtering step. For gathering, you can specify genes by
#' name or supply CES-style variant IDs to get all variants matching either. Or, leave
#' both empty and all available variants will be gathered. For filtering, use the granges
#' parameter to filter out variants that don't intersect an input GRanges object. Use
#' min_freq to filter out variants by frequency in the MAF data. 
#' 
#' 
#' Note that while intergenic SNVs have their nearest genes annotated in the SNV tables,
#' these variants will not be captured by gene-based selection with this function, since
#' they're not actually in any gene.
#' 
#' Definitions of some less self-explanatory columns:
#' \itemize{
#'   \item variant_name: short, often but not necessarily uniquely identifying name (use
#'   variant_id to guarantee uniqueness) 
#'   \item start/end: lowest/highest genomic positions overlapping variant
#'   \item variant_id: unique IDs for variants given the associated genome assembly version and the transcript data
#'   \item ref/alt: genomic reference and alternate alleles (for genomic
#'   variants; NA for AACs) 
#'   \item gene: the affected gene in AACs; for SNVs, the
#'   overlapping gene (or an abitrary gene if more than one overlaps), or the nearest gene
#'   for intergenic SNVs
#'   \item strand: for AACs, 1 if the reference sequence strand is the coding strand; -1 otherwise
#'   \item essential_splice: Variant is 1,2 bp upstream or
#'   1,2,5 bp downstream of an annotated splice position (edge case: if an SNV has
#'   multiple gene/transcript annotations, this doesn't say which one it's essential for)
#'   \item intergenic: variant does not overlap any gene ranges in the reference data
#'   \item trinuc_mut: for SNVs, the reference trinucleotide context, in deconstructSigs notation
#'   \item coding_seq: coding strand nucleotides in order of transcription
#'   \item center_nt_pos: regardless of strand, start/end give positions of two out of three AAC nucleotides; this
#'                       gives the position of the center nucleotide (maybe useful if the AAC spans a splice site)
#'   \item constituent_snvs: all SNVs that can produce a given variant 
#'   \item multi_anno_site: T/F whether variant has multiple gene/transcript/AAC annotations
#'   \item all_aac: For SNVs, all AACs that the variant is annotated as producing; for
#'   AACs, all of the AACs that are annotated for all of the constituent SNVs
#'   \item all_genes: all genes overlapping the variant in reference data
#'   \item covered_in: the names of all "covered regions" sets in the CESAnalysis that have coverage at the variant site
#'   \item maf_freq: number of occurrences of the variant in MAF data
#' }
#' 
#' @param cesa CESAnalysis with MAF data loaded and annotated (e.g., with \code{load_maf()})
#' @param genes include variants (coding and noncoding) within these genes
#' @param variant_ids include variants by ID: vectors of snv_id and/or aac_id, or short
#'   (potentially non-unique) names like KRAS_G12C (or "KRAS G12C")
#' @param granges filter out any variants not within input GRanges
#' @param min_freq filter out variants that occur in fewer than this many samples in the
#'   MAF data Defaults to 1; set to 0 to include all passing variants that are in the
#'   annotation tables. Note that variants that are not in the annotation tables will
#'   never be returned. Use \code{generate_variants()} to annotate variants absent from the
#'   data.
#' @param include_subvariants Some mutations "contain" other mutations. For example, in
#'   cancereffectsizeR's default hg19 reference data, KRAS_Q61H contains two constituent
#'   SNVs that both cause the same amino acid change: 12:25380275_T>G and 12:25380275_T>A.
#'   When include_subvariants = F (the default), and genes = "KRAS", output will be
#'   returned for KRAS_Q61H but not for the two SNVs (although their IDs will appear in
#'   the Q61H output). Set to true, and all three variants will be included in output, 
#'   assuming they don't get filtered out by other other options, like min_freq
#' @param gr_padding add +/- this many bases to every range specified in "granges"
#'   (stopping at chromsome ends, naturally)
#' @param ids_only instead of outputting a table of variant information, return a list
#' of selected variant IDs organized by variant type, without any annotations
#' @param collapse_lists Some output columns may have multiple elements per variant row.
#'   For example, all_genes may include multiple genes. When collapse_lists = T, these
#'   columns have their values collapsed into single comma-delimited strings (e.g.
#'   "Gene1,Gene2"), which are good for visual inspection of data. When collapse_lists =
#'   F, these columns are returned as "list columns", each element of which is a
#'   variable-length vector. These are better for advanced scripting, but the syntax to
#'   manipulate them is a little tricky, so the default is to collapse.
#' @return a data table with info on selected variants (see details), or a list of IDs
#' @export
select_variants = function(cesa, genes = NULL, variant_ids = NULL, granges = NULL, min_freq = 0, 
                           include_subvariants = F, gr_padding = 0, ids_only = F, collapse_lists = T) {
  
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis object")
  }
  bsg = get_cesa_bsg(cesa)
  
  if(cesa@maf[, .N] == 0) {
    stop("There are no variants in the CESAnalysis.")
  }
  if(length(cesa@mutations) == 0) {
    stop("There are no variant annotations (run annotate_variants or re-run load_maf with annotate = TRUE).")
  }
  if(! is.numeric(min_freq) || length(min_freq) > 1 || min_freq - as.integer(min_freq) != 0) {
    stop("min_freq should be 1-length integer at least zero")
  }
  if(! is.null(genes) & ! is(genes, "character")) {
    stop("genes should be character vector of genes to include or NULL (for no gene filtering)")
  }
  
  if(! is.null(variant_ids) & ! is(variant_ids, "character")) {
    stop("variant_ids should be character vector of variant_ids to include or NULL (for no variant ID filtering)")
  }
  
  if(! is.logical(include_subvariants) | length(include_subvariants) != 1) {
    stop("include_subvariants should be T/F")
  }
  
  if (! is.logical(collapse_lists) || length(collapse_lists) != 1) {
    stop("collapse_lists should be T/F")
  }
  
  if (! is.logical(ids_only) || length(ids_only) != 1) {
    stop("ids_only should be T/F")
  }
  
  if (ids_only == TRUE & collapse_lists == FALSE) {
    stop("There is no use setting collapse_list = FALSE when ids_only = TRUE")
  }
  
  
  selected_snv_ids = character()
  selected_aac_ids = character()
  
  if(! is.null(granges)) {
    if(! is(granges, "GRanges")) {
      stop("granges should be a GRanges object", call. = F)
    }
    granges = clean_granges_for_bsg(bsg = bsg, gr = granges, padding = gr_padding)
    genome_info = GenomeInfoDb::seqinfo(bsg)
    mutations_gr = GenomicRanges::makeGRangesFromDataFrame(cesa@mutations$snv, seqnames.field = "chr", start.field = "pos", 
                                                         end.field = "pos", seqinfo = genome_info)
    captured = cesa@mutations$snv[IRanges::overlapsAny(mutations_gr, granges, type = "within")]
    if (captured[, .N] == 0) {
      stop("No mutations captured by input granges.", call. = F)
    }
    gr_passing_snv_id = captured$snv_id
    gr_passing_aac_id = captured[!is.na(assoc_aac), unique(unlist(assoc_aac))]
  }
  
  if (length(variant_ids) > 0) {
    variant_ids = unique(variant_ids)
    matching_snv_ids = cesa@mutations$snv[variant_ids, snv_id, nomatch = NULL]
    matching_aac_ids = cesa@mutations$amino_acid_change[variant_ids, aac_id, nomatch = NULL]
    missing_ids = setdiff(variant_ids, c(matching_snv_ids, matching_aac_ids))
    
    # if any IDs are missing, try to interpret them as "short" AAC names (i.e., without protein ID)
    if (length(missing_ids) > 0) {
      missing_ids = gsub(' ', '_', missing_ids)
      tmp = cesa@mutations$amino_acid_change[, .(aac_id, variant_name = paste(gene, aachange, sep = "_"))]
      aac_matches = tmp[missing_ids, on = "variant_name"]
      
      missing_ids = aac_matches[is.na(aac_id), variant_name]
      if (length(missing_ids) > 0) {
        num_missing = length(missing_ids)
        missing_ids = paste(missing_ids, collapse = ", ")
        stop(num_missing, " variants given in variant_ids couldn't be found. Either they're bad IDs or there are no annotations for them.\n",
             missing_ids)
      }
      matching_aac_ids = c(matching_aac_ids, aac_matches$aac_id)
      if (any(duplicated(aac_matches$variant_name))) {
        msg = paste0("Shorthand amino-acid-change names (styled like \"KRAS_G12C\") were recognized and matched, ",
                     "with cancereffectsizeR-style aac_ids. However, some of your variant names matched more than ",
                     "one aac_id (e.g., same positional change on multiple protein isoforms), and all of these are included in output.")
        pretty_message(msg, black = F)
      } else {
        msg = paste0("Shorthand amino-acid-change names (styled like \"KRAS_G12C\") were recognized and uniquely ",
                     "paired with cancereffectsizeR's aac_ids.")
        pretty_message(msg)
      }
    }
    selected_with_variant_id = c(matching_snv_ids, matching_aac_ids) # save this for later in case there's a gr filter
    selected_snv_ids = matching_snv_ids
    selected_aac_ids = matching_aac_ids
  }
  
  
  if (length(genes) > 0) {
    genes = unique(genes)
    gene_names = get_ref_data(cesa, "gene_names")
    invalid_genes = genes[! genes %in% gene_names]
    num_invalid_genes = length(invalid_genes)
    if (num_invalid_genes > 0) {
      stop("Some of the selected genes do not appear in the CESAnalysis reference data:\n", paste(invalid_genes, collapse = ", "))
    }
    aac_in_genes = cesa@mutations$amino_acid_change[gene %in% genes, aac_id]
    
    # Note we're not returning intergenic SNVs that just have one of the chosen gene as their nearest gene
    genes_by_snv = cesa@mutations$snv[intergenic == FALSE, .(gene = unlist(genes)), by = "snv_id"]
    snv_in_genes = genes_by_snv[gene %in% genes, unique(snv_id)]
    
    if (include_subvariants == FALSE) {
      subvariant_snv_ids = cesa@mutations$amino_acid_change[aac_in_genes, unlist(constituent_snvs)]
      snv_in_genes = setdiff(snv_in_genes, subvariant_snv_ids)
    }
    selected_snv_ids = union(selected_snv_ids, snv_in_genes)
    selected_aac_ids = union(selected_aac_ids, aac_in_genes)
  }
  
  if (is.null(variant_ids) & is.null(genes)) {
    selected_aac_ids = cesa@mutations$amino_acid_change$aac_id
    selected_snv_ids = cesa@mutations$snv$snv_id
    if (include_subvariants == FALSE) {
      subvariant_snv_ids = cesa@mutations$amino_acid_change[selected_aac_ids, unlist(constituent_snvs)]
      selected_snv_ids = setdiff(selected_snv_ids, subvariant_snv_ids)
    }
  }
  
  if (! is.null(granges)) {
    selected_snv_ids = intersect(selected_snv_ids, gr_passing_snv_id)
    selected_aac_ids = intersect(selected_aac_ids, gr_passing_aac_id) 
    if (! is.null(variant_ids) & any(! selected_with_variant_id %in% c(selected_snv_ids, selected_aac_ids))) {
      message("FYI, some of the variants specified with variant_ids were filtered out by your input GRanges.")
    }
  }
  
  selected_snv = cesa@mutations$snv[selected_snv_ids]
  setkey(selected_snv, "snv_id")
  selected_aac = cesa@mutations$amino_acid_change[selected_aac_ids]
  setkey(selected_aac, "aac_id")
  
  # Tabulate variants in MAF data and apply frequency filter
  snv_counts = cesa@maf[variant_type == "snv", .N, by = "variant_id"][N >= min_freq]
  aac_counts = cesa@maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac))][, .N, by = "aac_id"][N >= min_freq]
  selected_snv[, maf_frequency := 0]
  selected_snv = selected_snv[snv_counts, maf_frequency := N][maf_frequency >= min_freq]
  selected_aac[, maf_frequency := 0]
  selected_aac = selected_aac[aac_counts, maf_frequency := N][maf_frequency >= min_freq]
  
  # Return just the IDs if requested
  if (ids_only) {
    return(list(aac = selected_aac$aac_id, snv = selected_snv$snv_id))
  }
  
  # Annotate SNV table and prepare to merge with AACs
  selected_snv[, variant_type := "snv"]
  selected_snv[, variant_name := snv_id] # SNV IDs are already short and uniquely identifying
  selected_snv[, constituent_snvs := list(NA_character_)]
  selected_snv[, strand := NA_integer_] # because AAC table is +1/-1
  selected_snv[, c("start", "end") := .(pos, pos)]
  selected_snv[, pos := NULL]
  # for convenience, take 1 gene per SNV for the output gene column (all will appear in all_genes column)
  if (selected_snv[, .N] > 0) {
    selected_snv[, gene := sapply(genes, function(x) x[1])] 
  } else {
    selected_snv[, gene := character()]
  }
  selected_snv[intergenic == T, genes := list(NA_character_)]
  selected_snv[, multi_anno_site := F][which(sapply(assoc_aac, length) > 1 | sapply(genes, length) > 1), multi_anno_site := T]
  setnames(selected_snv, c("genes", "assoc_aac", "snv_id"), c("all_genes", "all_aac", "variant_id"))
  # AACs get a short variant name that might not be uniquely identifying if a gene has more than one CDS
  selected_aac[, variant_name := paste(gene, aachange, sep = "_")]
  selected_aac[, variant_type := "aac"]
  selected_aac[, intergenic := FALSE]
  selected_aac[, start := pmin(nt1_pos, nt3_pos)]
  selected_aac[, end := pmax(nt1_pos, nt3_pos)]
  selected_aac[, center_nt_pos := nt2_pos]
  selected_aac[, c("nt1_pos", "nt2_pos", "nt3_pos") := NULL]
  
  if (selected_aac[, .N] > 0) {
    aac_to_snv = selected_aac[, .(snv_id = unlist(constituent_snvs)), by = "aac_id"]
    aac_to_snv[, c("genes", "assoc_aac") := cesa@mutations$snv[aac_to_snv$snv_id, .(genes, assoc_aac)]]
    all_genes_by_aac_id = aac_to_snv[, .(genes = .(unique(unlist(genes)))), by = "aac_id"]
    assoc_aac_by_aac_id = aac_to_snv[, .(assoc_aac = .(unique(unlist(assoc_aac)))), by = "aac_id"]
    selected_aac[all_genes_by_aac_id, all_genes := genes, on = "aac_id"]
    selected_aac[assoc_aac_by_aac_id, all_aac := assoc_aac, on = "aac_id"]
  } else {
    selected_aac[, all_genes := list(NA_character_)]
    selected_aac[, all_aac := list(NA_character_)]
  }
  selected_aac[, multi_anno_site := F][which(sapply(all_aac, length) > 1 | sapply(all_genes, length) > 1), multi_anno_site := T]
  setnames(selected_aac, "aac_id", "variant_id")

  # Combine SNV and AAC tables
  combined = rbindlist(list(selected_aac, selected_snv), use.names = T, fill = T)
  if(combined[, .N] == 0) {
    message("No variants passed selection criteria!")
    return(NULL)
  }
  
  all_cov_cols = character() # for output column name ordering
  for (curr_group in cesa@groups) {
    curr_samples = cesa@samples[group == curr_group]
    num_wgs_samples = curr_samples[covered_regions == "genome", .N]
    cov_counts = curr_samples[, .N, by = "covered_regions"]
    
    # this should be made more elegant at some point
    unique_combos = setdiff(unique(combined$covered_in), c(list(character()), list(NULL))) # various empty entries in sites just covered in WGS
    setkey(cov_counts, "covered_regions")
    
    # If all full-coverage WGS data, then there will be no unique_combos
    cov_count_colname = ifelse(length(cesa@groups) == 1, "samples_covering", paste0("samples_covering_in_", curr_group))
    all_cov_cols = c(all_cov_cols, cov_count_colname)
    if (length(unique_combos) > 0) {
      combo_counts = sapply(unique_combos, function(x) sum(cov_counts[x, N], na.rm = T))
      names(combo_counts) = sapply(unique_combos, function(x) paste0(x, collapse = "_"))
      combined[, (cov_count_colname) := combo_counts[S4Vectors::unstrsplit(covered_in, sep = "_")]]
      repl = combined[[cov_count_colname]]
      repl[is.na(repl)] = 0
      combined[, (cov_count_colname) := repl]
      combined[, (cov_count_colname) := combined[[cov_count_colname]] + num_wgs_samples]
    } else {
      combined[, (cov_count_colname) := num_wgs_samples]
    }
  }
  
  if (length(cesa@groups) > 1) {
    combined[, total_samples_covering := rowSums(.SD), .SDcols = all_cov_cols]
    all_cov_cols = c(all_cov_cols, "total_samples_covering")
  }
  # convert 0-length covered_in to NA
  combined[which(sapply(covered_in, length) == 0), covered_in := list(NA_character_)]
  
  
  # Break down frequency counts
  maf_freq_cols = character()
  if (length(cesa@groups) > 1) {
    for (curr_group in cesa@groups) {
      curr_col = paste0("maf_freq_in_", curr_group)
      maf_freq_cols = c(maf_freq_cols, curr_col)
      group_maf = cesa@maf[Unique_Patient_Identifier %in% cesa@samples[group == curr_group, Unique_Patient_Identifier]]
      combined[, (curr_col) := 0]
      snv_counts = group_maf[variant_type == "snv", .N, by = "variant_id"][N >= min_freq]
      if(snv_counts[, .N] > 0) {
        combined[snv_counts, (curr_col) := N, on = "variant_id"]
        # can't be AACs unless there are SNVs, hence nested
        aac_counts = group_maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac)), by = "variant_id"]
        if(aac_counts[, .N] > 0) {
          aac_counts = aac_counts[, .N, by = "aac_id"][N >= min_freq]
          combined[aac_counts, (curr_col) := N, on = c(variant_id = "aac_id")]
        }
      }
    }
    setnames(combined, "maf_frequency", "total_maf_freq")
    maf_freq_cols = c(maf_freq_cols, "total_maf_freq")
  } else {
    maf_freq_cols = "maf_frequency"
  }
  # collapse list columns, if specified
  if (collapse_lists) {
    # Problem: unstrsplit converts NA to "NA"
    # "NA" is not a valid value for any of these except all_genes, and going to assume there is not a gene called "NA"
    list_cols = c("constituent_snvs", "covered_in", "all_genes", "all_aac")
    combined[, (list_cols) := lapply(.SD, function(x) S4Vectors::unstrsplit(x, sep = ",")), .SDcols = list_cols]
    combined[, (list_cols) := lapply(.SD, function(x) gsub('NA', NA_character_, x)), .SDcols = list_cols]

  }

  # order output in chr/pos order
  setkey(combined, "chr")
  combined = combined[BSgenome::seqnames(bsg), nomatch = NULL][, .SD[order(start)], by = "chr"]
  setcolorder(combined, c("variant_name", "variant_type", "chr", "start", "end", "variant_id", "ref", "alt", "gene", 
                          "strand", "aachange", "essential_splice", "intergenic", "trinuc_mut", "aa_ref", "aa_pos", "aa_alt", "coding_seq", 
                          "center_nt_pos", "pid", "constituent_snvs", "multi_anno_site", "all_aac", "all_genes",
                          "covered_in", unlist(zipup(all_cov_cols, maf_freq_cols))))
  
  # notify about sites with multiple genes/CDS
  if(any(combined$multi_anno_site)) {
    gene_note = ifelse(is.null(genes), "", ", regardless of your \"genes\" input")
    collapse_note = ifelse(collapse_lists == T,  'comma-delimited', 'listed')
    msg = paste0("Note: Some of your variants cover sites that overlap more than one gene, transcript or coding sequence. These are ",
                 collapse_note, " in the \"all_genes\" and \"all_aac\" columns. (For SNVs, the gene column gives a single gene arbitrarily", 
                 gene_note, ".) You can find these variants by filtering on the multi_anno_site column.")
    pretty_message(msg)
  }
  return(combined[]) # brackets force the output to print when unassigned (should automatically, but this is a known data.table issue)
}


#' add_covered_regions
#' @param target_cesa CESAnalysis with annotated variants that the covered regions will be added to
#' @param source_cesa Another CESAnalysis to copy all covered regions from
#' @param covered_regions A GRanges object or BED file path with genome build matching the target_cesa,
#'                      if not using source_cesa
#' @param covered_regions_name A name to identify the covered regions, if not using source_cesa
#' @param coverage_type exome or targeted, if not using source_cesa
#' @param covered_regions_padding optionally, add +/- this many bp to each interval in covered_regions
#' @return CESAnalysis given in target_cesa, with the new covered regions added
#' @export
add_covered_regions = function(target_cesa = NULL, source_cesa = NULL, covered_regions = NULL, 
                               covered_regions_name = NULL, coverage_type = NULL, covered_regions_padding = 0) {
  prev_recording_status = target_cesa@advanced$recording
  target_cesa = update_cesa_history(target_cesa, match.call())
  target_cesa@advanced$recording = F
  if (! is(target_cesa, "CESAnalysis")) {
    stop("target_cesa should be a CESAnalysis", call. = F)
  }
  if (! is(covered_regions_padding, "numeric") || length(covered_regions_padding) > 1 || covered_regions_padding < 0 ||
      covered_regions_padding - as.integer(covered_regions_padding) != 0) {
    stop("covered_regions_padding should be 1-length integer", call. = F)
  }
  
  coverage_args = list(covered_regions, covered_regions_name, coverage_type)
  # Must supply just source_cesa or covered regions information
  if (is.null(source_cesa) & any(sapply(coverage_args, is.null))) {
    stop("To add a covered regions set, you need to supply covered_regions, covered_regions_name, and coverage_type.")
  }
  if(! is.null(source_cesa) & ! all(sapply(coverage_args, is.null))) {
    stop("Use source_cesa to copy all covered regions from another CESAnalysis, or\n",
         "covered_regions, covered_regions_name, and coverage_type to add one covered regions set.")
  }
  
  # Validate both possible sets of arguments
  if (! is.null(source_cesa)) {
    if (! is(source_cesa, "CESAnalysis")) {
      stop("source_cesa should be a CESAnalysis")
    }
    
    if (! identical(target_cesa@advanced$annotated, T)) {
      stop("target_cesa should be annotated")
    }
    if(! identical(get_cesa_bsg(target_cesa), get_cesa_bsg(source_cesa))) {
      stop("target_cesa and source_cesa have different reference genomes")
    }
    if (covered_regions_padding != 0) {
      stop("covered_regions_padding can't be used with source_cesa (leave it unspecified).")
    }
    exome_sets = names(source_cesa@coverage$exome)
    if (any(c("exome", "exome+") %in% exome_sets)) {
      warning("Generic exome/exome+ covered regions sets found in source_cesa. These will not be copied to the target\n",
              "CESAnalysis. If you really want to copy them over, extract their GRanges from the CESAnalysis and add\n",
              "them individually. (As the \"exome\" set is a set of default coverage intervals for generic exome\n",
              "data, and \"exome+\" is an expanded version of that set, specific to a CESAnalysis, that covers all\n",
              "variants loaded as generic exome data, they're probably not useful to copy.)")
      exome_sets = setdiff(exome_sets, c("exome", "exome+"))
    }
    for (exome_set in exome_sets) {
      target_cesa = assign_gr_to_coverage(target_cesa, gr = source_cesa@coverage$exome[[exome_set]], 
                                          covered_regions_name = exome_set, coverage_type = "exome")
    }
    tgs_sets = names(source_cesa@coverage$targeted)
    for (tgs_set in tgs_sets) {
      target_cesa = assign_gr_to_coverage(target_cesa, gr = source_cesa@coverage$targeted[[tgs_set]],
                                          covered_regions_name = tgs_set, coverage_type ="targeted")
    }
    target_cesa@advanced$recording = prev_recording_status
    return(update_covered_in(target_cesa))
  } else {
    target_cesa@advanced$recording = prev_recording_status
    return(.add_covered_regions(cesa = target_cesa, coverage_type = coverage_type, covered_regions = covered_regions,
                     covered_regions_name = covered_regions_name,covered_regions_padding = covered_regions_padding, 
                     update_anno = target_cesa@advanced$annotated))
  }
}

#' .add_covered_regions
#' @param covered_regions A GRanges object or BED file path with genome build matching the target_cesa,
#'                      if not using source_cesa
#' @param covered_regions_name A name to identify the covered regions, if not using source_cesa
#' @param coverage_type exome or targeted, if not using source_cesa
#' @param covered_regions_padding optionally, add +/- this many bp to each interval in covered_regions
#' @param update_anno T/F, whether to update the covered_in fields in variant annotations
#' @keywords internal
#' @return CESAnalysis given in target_cesa, with the new covered regions added
.add_covered_regions = function(cesa, coverage_type, covered_regions, covered_regions_name, covered_regions_padding, update_anno) {
  if (! is.character(coverage_type) | length(coverage_type) != 1 | ! coverage_type %in% c("exome", "targeted")) {
    stop("coverage_type should be exome or targeted.", call. = F)
  }
  if (! is.character(covered_regions_name) | length(covered_regions_name) != 1 | is.na(covered_regions_name)) {
    stop("covered_regions_name should be 1-length character.")
  }
  # exome, exome+, genome are reserved, and we don't want "NA" to be the name
  if (tolower(covered_regions_name) %in% c("exome", "exome+", "genome", "na")) {
    stop("Please pick a different covered_regions_name; you chose one reserved for internal use.")
  }
  # covered_regions_name must start with letter, contain only letters, numbers, underscores, hyphen, period
  legal_name = '^[a-z][0-9a-z\\_\\-\\.]*$'
  if (! grepl(legal_name, tolower(covered_regions_name), perl = T)) {
    stop("Invalid covered_regions_name. The name must start with a letter and contain only letters, numbers, and '-', '_', '.'.")
  }
  bad_covered_regions_msg = "covered_regions should be GRanges or BED file path."
  if (is.character(covered_regions)) {
    if (length(covered_regions) != 1) {
      stop(bad_covered_regions_msg, call. = F)
    }
    if (! file.exists(covered_regions)) {
      stop("BED file not found; check path?", call. = F)
    }
    gr = rtracklayer::import.bed(covered_regions)
  } else if (is(covered_regions, 'GRanges')) {
    gr = covered_regions
  } else {
    stop(bad_covered_regions_msg, call. = F)
  }
  gr = clean_granges_for_bsg(bsg = get_cesa_bsg(cesa), gr = gr, padding = covered_regions_padding)
  cesa = assign_gr_to_coverage(cesa, gr = gr, covered_regions_name = covered_regions_name, coverage_type = coverage_type)
  if (update_anno) {
    cesa = update_covered_in(cesa)
  }
  return(cesa)
}


#' assign_gr_to_coverage
#' 
#' Adds a validated GRanges object as a CESAnalysis's coverage set. Called by 
#' add_covered_regions() after various checks pass.
#' 
#' Special handling occurs if covered_regions_name is "exome+".
#' @param cesa CESAnalysis to receive the gr
#' @param gr GRanges
#' @param covered_regions_name unique name for the covered regions
#' @param coverage_type "exome" or "targeted"
#' @keywords internal
assign_gr_to_coverage = function(cesa, gr, covered_regions_name, coverage_type) {
  # if covered_regions_name was already used in a previous load_maf call, it must have
  # have been the same data type (exome, targeted), and the grs must match exactly
  other_coverage_types = setdiff(c("exome", "targeted", "genome"), coverage_type)
  if (covered_regions_name %in% unlist(lapply(cesa@coverage[other_coverage_types], names))) {
    stop("The covered_regions_name (", covered_regions_name, ")", " has already been used for a different type of sequencing data.")
  } 
  else if (covered_regions_name %in% names(cesa@coverage[[coverage_type]])) {
    if (! identical(gr, cesa@coverage[[coverage_type]][[covered_regions_name]])) {
      stop("MAF data was previously loaded in using the same covered_regions_name (", covered_regions_name, "),\n",
           "but the covered_regions do not exactly match. Perhaps the input BED files (or GRanges) are from\n",
           "different sources, or different amounts of interval padding were used.")
    } else {
      # nothing more to do if the covered regions are already present
      return(cesa)
    }
  }
  
  # If possible, see if covered regions size resembles exome data
  if (coverage_type == "exome" & check_for_ref_data(cesa, "generic_exome_gr")) {
    covered_regions_bases_covered = sum(IRanges::width(IRanges::ranges(gr)))
    generic_bases_covered = sum(IRanges::width(IRanges::ranges(get_ref_data(cesa, "generic_exome_gr"))))
    if (covered_regions_bases_covered / generic_bases_covered < .4) {
      warning(paste0("Input coverage ranges are described as exome but are less than 40% of the size of this genome's default exome intervals.\n",
                     "This might make sense if your exome capture array is  lean, but if this is actually targeted sequencing data,\n",
                     "start over with the coverage=\"targeted\"."))
    }
  }
  
  if (! coverage_type %in% names(cesa@coverage)) {
    cesa@coverage[[coverage_type]] = list()
  }
  cesa@coverage[[coverage_type]][[covered_regions_name]] = gr
  return(cesa)
}



#' Add variant annotations
#' 
#' Use this function to add variant annotations to your CESAnalysis, either using
#' a vector of CES-style snv_ids or a source CESAnalysis. Either way, the variants 
#' will be re-annotated using the target CESAnalysis's associated reference data. In other words,
#' when you use source_cesa, only the snv_ids are extracted, not any other information, in order
#' to avoid ugly situations involving conflicting reference data within a CESAnalysis. 
#' @param target_cesa CESAnalysis to receive variant annotations
#' @param snv_id character vector of CES-style snv_ids to validate, annotate (including
#'   finding and adding associated amino-acid-change mutations), and add to the target
#'   CESAnalysis
#' @param source_cesa source CESAnalysis whose snv_ids will be used to create new
#'   annotations in the target CESAnalysis
#' @export
add_variants = function(target_cesa = NULL, snv_id = NULL, source_cesa = NULL) {
  if(! is(target_cesa, "CESAnalysis")) {
    stop("target_cesa should be a CESAnalysis", call. = F)
  }
  target_cesa = update_cesa_history(target_cesa, match.call())
  if(! identical(target_cesa@advanced$annotated, T)) {
    stop("target_cesa should be annotated", call. = F)
  }

  if(! xor(is.null(snv_id), is.null(source_cesa))) {
    stop("Add annotations via snv_id or source_cesa, but not both", call. = F)
  }

  if (! is.null(source_cesa)) {
    if(! is(source_cesa, "CESAnalysis")) {
      stop("source_cesa should be a CESAnalysis", call. = F)
    }
    source_snv_table = source_cesa@mutations$snv
    if (is.null(source_snv_table)) {
      stop("source_cesa has no SNV annotations", call. = F)
    }
    snvs_to_annotate = source_cesa@mutations$snv$snv_id
  }
  
  # Load reference data if not already present
  if (! target_cesa@ref_key %in% ls(.ces_ref_data)) {
    preload_ref_data(target_cesa@ref_data_dir)
  }
  
  if(! is.null(snv_id)) {
    if(! is(snv_id, "character") | length(snv_id) == 0) {
      stop("Expected snv_id to be character vector of snv_ids (e.g., 1:100_A>G", call. = F)
    }
    # will stop with errror if any IDs fail validation
    validate_snv_ids(snv_id, .ces_ref_data[[target_cesa@ref_key]]$genome)
    snvs_to_annotate = snv_id
  }

  maf = as.data.table(tstrsplit(snvs_to_annotate, split = '[:_>]'))
  colnames(maf) = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
  maf$Unique_Patient_Identifier = NA_character_
  setcolorder(maf, "Unique_Patient_Identifier")
  maf$variant_type = "snv"
  maf[, Start_Position := as.numeric(Start_Position)]
  
  target_cesa@maf = rbind(target_cesa@maf, maf, fill = T)
  prev_recording_status = target_cesa@advanced$recording
  target_cesa@advanced$recording = F
  target_cesa = annotate_variants(target_cesa)
  target_cesa@advanced$recording = prev_recording_status
  target_cesa@maf = target_cesa@maf[! is.na(Unique_Patient_Identifier)]
  return(target_cesa)
}

#' validate_snv_ids
#' 
#' Ensures SNV IDs are valid for the given genome
#' 
#' @param snv_ids character vector of snv_ids
#' @param bsg BSgenome for getting reference sequence
#' @keywords internal
validate_snv_ids = function(snv_ids, bsg) {
  dt = as.data.table(tstrsplit(snv_ids, split = '[:_>]'))
  nt = c("A", "C", "G", "T")
  if(length(dt) != 4 | anyNA(dt)) {
    stop("One of the input IDs does not match snv_id format (e.g., 1:100_A>G)")
  }
  names(dt) = c("chr", "start", "ref", "alt")
  
  fail_msg = "Check that chromosome names and genome assembly match the CESAnalysis.\nOriginal error/warning:"
  tryCatch(
    {
      seqs = as.character(getSeq(bsg, makeGRangesFromDataFrame(dt, start.field = "start", end.field = "start")))
    },
    error = function(e) {
      stop(paste(fail_msg, e, sep = "\n"), call. = F)
    },
    warning = function(e) {
      stop(paste(fail_msg, e, sep = "\n"), call. = F)
    }
  )
  
  if(! all(seqs == dt$ref)) {
    stop("Incorrect reference allele in one or more SNV IDs.")
  }
  
  if(! all(dt$alt %in% nt)) {
    stop("SNV alt alleles are not all single DNA bases")
  }
  if (! all(dt$alt != dt$ref)) {
    stop("Some SNV alt alleles match the ref alleles.")
  }
}

