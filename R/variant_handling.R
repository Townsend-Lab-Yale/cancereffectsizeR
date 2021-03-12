#' Select and filter variants
#'
#' This function helps you find and view variant data from your CESAnalysis's MAF data and
#' mutation annotation tables. By default, almost all amino-acid-change mutations and
#' noncoding SNVs are returned. You can apply a series of filters to restrict output to
#' certain genes or genomic regions or require a minimum variant frequency in MAF data.
#' You can also specify some variants to include in output regardless of filters with
#' \code{variant_passlist}. Special behavior: If \code{variant_passlist} is used by
#' itself, then only those specified variants will be returned.
#' 
#' Only variants that are present in the CESAnalysis's annotation tables can be returned,
#' which by default are those present in the MAF data. To select variants absent from MAF
#' data, you must first call add_variants() to add them to the CESAnalysis. Note that
#' while intergenic SNVs have their nearest genes annotated in the SNV tables, these
#' variants will not be captured by gene-based selection with this function, since they're
#' not actually in any gene.
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
#'   overlapping gene (or an arbitrary gene if more than one overlaps), or the nearest gene
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
#' @param genes Filter variants to specified genes.
#' @param min_freq Filter out variants with MAF frequency below threshold (default 0).
#'   Note that variants that are not in the annotation tables will never be returned. Use
#'   \code{add_variants()} to include variants absent from MAF data in your CESAnalysis.
#' @param variant_passlist Vector of variant IDs to include in output regardless of
#'   filtering options. You can use CES-style AAC and SNV IDs or variant names like 
#'   "KRAS G12C". If this argument is used by itself (without any filtering arguments),
#'   then only these specified variants will be returned.
#' @param gr Filter out any variants not within input GRanges +/- \code{padding} bases.
#' @param variant_position_table Filter out any variants that don't intersect the
#'   positions given in chr/start/end of this table (1-based closed coordinates).
#'   Typically, the table comes from a previous \code{select_variants} call and can be expanded
#'   with \code{padding}. (Gritty detail: Amino acid change SNVS get special handling. Only the
#'   precise positions in start, end, and center_nt_pos are used. This avoids intersecting
#'   extra variants between start/end, which on splice-site-spanning variants can be many
#'   thousands.)
#' @param padding add +/- this many bases to every range specified in \code{gr} or
#'   \code{variant_position_table} (stopping at chromosome ends, naturally).
#' @param collapse_lists Some output columns may have multiple elements per variant row.
#'   For example, all_genes may include multiple genes. These variable-length vectors
#'   allow advanced filtering and manipulation, but the syntax can be tricky. Optionally,
#'   set collapse_lists = T to convert these columns to comma-delimited strings, which are
#'   sometimes easier to work with.
#' @param include_subvariants Some mutations "contain" other mutations. For example, in
#'   cancereffectsizeR's ces.refset.hg19, KRAS_Q61H contains two constituent
#'   SNVs that both cause the same amino acid change: 12:25380275_T>G and 12:25380275_T>A.
#'   When include_subvariants = F (the default), and genes = "KRAS", output will be
#'   returned for KRAS_Q61H but not for the two SNVs (although their IDs will appear in
#'   the Q61H output). Set to true, and all three variants will be included in output, 
#'   assuming they don't get filtered out by other other options, like min_freq. If you
#'   set this to TRUE, you can't directly plug the output table into selection functions.
#'   However, you can pick a non-overlapping set of variant IDs from the output table 
#'   and re-run \code{select_variants()} to put those variants into a new table for
#'   selection functions.
#' @param remove_secondary_aac Default TRUE, except overridden (effectively FALSE) when
#'   include_subvariants = T. Occasionally, due to overlapping coding sequence definitions
#'   in reference data, a site can have more than one amino-acid-change mutation
#'   annotation. To avoid returning the same genome-positional variants multiple times,
#'   the default is to return one AAC in these situations. Tiebreakers are MAF frequency,
#'   essential splice site status, non-silent status, gene mutation frequency,
#'   alphabetical. The omitted AACs will still be included in the all_aac column, and if
#'   any constituent SNVS of the AACs don't overlap, the non-overlapping SNVs from the
#'   excluded AACs will still be included in output. If you set remove_secondary_aac to
#'   FALSE, you can't put the output variant table in selection calculation functions. An
#'   alternative is to set to FALSE, pick which (non-overlapping) variants you want, and
#'   then re-run select_variants() with those variants specified in \code{variant_passlist}.
#' @return A data table with info on selected variants (see details), or a list of IDs.
#' @export
select_variants = function(cesa, genes = NULL, min_freq = 0, variant_passlist = NULL, gr = NULL, variant_position_table = NULL, 
                            include_subvariants = F, padding = 0, collapse_lists = F, remove_secondary_aac = TRUE) {
  
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis object")
  }
  bsg = get_cesa_bsg(cesa)
  
  if(length(cesa@mutations) == 0) {
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
  
  if(! is.null(variant_passlist) & ! is(variant_passlist, "character")) {
    stop("variant_passlist should be character vector of variants to include or NULL.")
  }
  
  if(! is.logical(include_subvariants) | length(include_subvariants) != 1) {
    stop("include_subvariants should be T/F")
  }
  
  if (! is.logical(collapse_lists) || length(collapse_lists) != 1) {
    stop("collapse_lists should be T/F")
  }
  
  if (! is.logical(remove_secondary_aac) || length(remove_secondary_aac) != 1) {
    stop("collapse_lists should be T/F")
  }
  if (include_subvariants) {
    remove_secondary_aac = FALSE
  }

  # collect all variants, unless just variant_passlist specified
  if (is.null(gr) && is.null(variant_position_table) && is.null(genes) && min_freq == 0 && length(variant_passlist) > 0) {
    selected_aac_ids = character()
    selected_snv_ids = character()
  } else {
    selected_aac_ids = cesa@mutations$amino_acid_change$aac_id
    selected_snv_ids = cesa@mutations$snv$snv_id
  }
  
  
  # handle variant_position_table or gr (for simplicity, not allowing both)
  final_gr = NULL
  if(! is.null(variant_position_table)) {
    if (! is.null(gr)) {
      stop("Sorry, you can't use both variant_position_table and gr. You should be able to get what you want by ",
           "running select_variants() sequentially.")
    }
    # get_gr_from_table handles validation of variant_position_table
    final_gr = clean_granges_for_cesa(cesa = cesa, gr = get_gr_from_table(variant_position_table), padding = padding)
  }
  
  if(! is.null(gr)) {
    if(! is(gr, "GRanges")) {
      stop("gr should be a GRanges object", call. = F)
    }
    gr = clean_granges_for_cesa(cesa = cesa, gr = gr, padding = padding)
    if (! is.null(final_gr)) {
      final_gr = GenomicRanges::intersect(final_gr, gr)
    } else {
      final_gr = gr
    }
  }
  
  # final_gr may derive from variant_position_table or gr
  if(! is.null(final_gr)) {
    genome_info = GenomeInfoDb::seqinfo(bsg)
    mutations_gr = GenomicRanges::makeGRangesFromDataFrame(cesa@mutations$snv, seqnames.field = "chr", start.field = "pos", 
                                                           end.field = "pos", seqinfo = genome_info)
    captured = cesa@mutations$snv[IRanges::overlapsAny(mutations_gr, final_gr, type = "within")]
    if (captured[, .N] == 0) {
      stop("No mutations captured by input genomic positions (gr/variant_position_table).", call. = F)
    }
    selected_snv_ids = intersect(selected_snv_ids, captured$snv_id)
    aac_passing_gr = union(character(), captured[!is.na(assoc_aac), unique(unlist(assoc_aac))])
    selected_aac_ids = intersect(selected_aac_ids, aac_passing_gr)
  }
  
  # Apply gene filter
  if (length(genes) > 0) {
    genes = unique(genes)
    gene_names = get_ref_data(cesa, "gene_names")
    invalid_genes = genes[! genes %in% gene_names]
    num_invalid_genes = length(invalid_genes)
    if (num_invalid_genes > 0) {
      stop("Some of the selected genes do not appear in the CESAnalysis reference data:\n", paste(invalid_genes, collapse = ", "))
    }
    
    aac_in_genes = cesa@mutations$amino_acid_change[gene %in% genes, aac_id]
    selected_aac_ids = intersect(selected_aac_ids, aac_in_genes)
    
    # Note we're not returning intergenic SNVs that just have one of the chosen genes as their nearest gene
    genes_by_snv = cesa@mutations$snv[intergenic == FALSE, .(gene = unlist(genes)), by = "snv_id"]
    snv_in_genes = genes_by_snv[gene %in% genes, unique(snv_id)]
    selected_snv_ids = intersect(selected_snv_ids, snv_in_genes)
  }
  
  # Include variants by ID
  passlisted_ids = character()
  if (length(variant_passlist) > 0) {
    variant_passlist = unique(variant_passlist)
    matching_snv_ids = cesa@mutations$snv[variant_passlist, snv_id, nomatch = NULL]
    matching_aac_ids = cesa@mutations$amino_acid_change[variant_passlist, aac_id, nomatch = NULL]
    missing_ids = setdiff(variant_passlist, c(matching_snv_ids, matching_aac_ids))
    
    # if any IDs are missing, try to interpret them as "short" AAC names (i.e., without protein ID)
    if (length(missing_ids) > 0) {
      missing_ids = gsub(' ', '_', missing_ids)
      tmp = cesa@mutations$amino_acid_change[, .(aac_id, variant_name = paste(gene, aachange, sep = "_"))]
      aac_matches = tmp[missing_ids, on = "variant_name"]
      
      missing_ids = aac_matches[is.na(aac_id), variant_name]
      if (length(missing_ids) > 0) {
        num_missing = length(missing_ids)
        missing_ids = paste(missing_ids, collapse = ", ")
        stop(num_missing, " variants given in variant_passlist couldn't be found. Either they're bad IDs or there are no annotations for them.\n",
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
    # under include_subvariants, all SNVs of passlisted AACs get included
    if (include_subvariants) {
      matching_snv_ids = union(matching_snv_ids, cesa@mutations$amino_acid_change[matching_aac_ids, unique(unlist(constituent_snvs))])
    }
    
    selected_snv_ids = union(selected_snv_ids, matching_snv_ids)
    selected_aac_ids = union(selected_aac_ids, matching_aac_ids)
    passlisted_ids = c(matching_snv_ids, matching_aac_ids)
  }
  
  # Remove variants contained in other variants
  # Rare exception: Constituent SNVs of AACs that get tossed in de-overlapping
  # need to be saved if they pass filters in their own right and don't overlap the AACs kep
  snvs_to_save = selected_snv_ids
  
  # We'll use the SNV counts to filter selected_snv_ids shortly)
  if(cesa@maf[, .N] > 0) {
    snv_counts = cesa@maf[variant_type == "snv", .N, by = "variant_id"][N >= min_freq]
    snvs_to_save = intersect(snvs_to_save, snv_counts$variant_id)
  } else {
    snv_counts = NA
    snvs_to_save = character()
  }

  if (include_subvariants == FALSE) {
    subvariant_snv_ids = cesa@mutations$amino_acid_change[selected_aac_ids, unlist(constituent_snvs), nomatch = NULL]
    selected_snv_ids = setdiff(selected_snv_ids, subvariant_snv_ids)
  }
  selected_snv = cesa@mutations$snv[selected_snv_ids]
  setkey(selected_snv, "snv_id")
  selected_aac = cesa@mutations$amino_acid_change[selected_aac_ids]
  setkey(selected_aac, "aac_id")

  # Tabulate variants in MAF data and apply frequency filter, exempting passlist variants from filters
  if (cesa@maf[, .N] > 0) {
    aac_counts = cesa@maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac))][, .N, by = "aac_id"][N >= min_freq]
  } else {
    aac_counts = NA
  }
  selected_snv[, maf_frequency := 0]
  selected_snv = selected_snv[snv_counts, maf_frequency := N]
  good_snv = union(passlisted_ids, selected_snv[maf_frequency >= min_freq, snv_id]) # okay to mix AAC/SNV in passlist
  selected_snv = selected_snv[good_snv, nomatch = NULL] # nomatch drops the AACs from passlist
  selected_aac[, maf_frequency := 0]
  selected_aac = selected_aac[aac_counts, maf_frequency := N]
  good_aac = union(passlisted_ids, selected_aac[maf_frequency >= min_freq, aac_id])
  selected_aac = selected_aac[good_aac, nomatch = NULL]
  if (any(selected_snv$maf_frequency < min_freq) || any(selected_aac$maf_frequency < min_freq)) {
    pretty_message("Note: Some of your passlist variants have MAF frequency < maf_freq. They will still appear in output.")
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
  if (cesa@samples[, .N] > 0) {
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
  }
  # convert 0-length covered_in to NA
  combined[which(sapply(covered_in, length) == 0), covered_in := list(NA_character_)]
  
  
  # Break down frequency counts
  maf_freq_cols = character()
  if (length(cesa@groups) > 1 & cesa@maf[, .N] > 0) {
    for (curr_group in cesa@groups) {
      curr_col = paste0("maf_freq_in_", curr_group)
      maf_freq_cols = c(maf_freq_cols, curr_col)
      group_maf = cesa@maf[Unique_Patient_Identifier %in% cesa@samples[group == curr_group, Unique_Patient_Identifier]]
      combined[, (curr_col) := 0]
      snv_counts = group_maf[variant_type == "snv", .N, by = "variant_id"]
      if(snv_counts[, .N] > 0) {
        combined[snv_counts, (curr_col) := N, on = "variant_id"]
        # can't be AACs unless there are SNVs, hence nested
        aac_counts = group_maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac)), by = "variant_id"]
        if(aac_counts[, .N] > 0) {
          aac_counts = aac_counts[, .N, by = "aac_id"]
          combined[aac_counts, (curr_col) := N, on = c(variant_id = "aac_id")]
        }
      }
    }
    setnames(combined, "maf_frequency", "total_maf_freq")
    maf_freq_cols = c(maf_freq_cols, "total_maf_freq")
  } else if(cesa@maf[, .N] > 0) {
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
  
  # handle overlapping  mutations
  if (remove_secondary_aac) {
    multi_hits = combined[sapply(all_aac, length) > 1 & variant_type == "aac"]
    num_to_check = multi_hits[, .N]
    if (num_to_check > 0) {
      setkey(multi_hits, "variant_id")
      # for tie-breaking, count how many mutations are in each gene found in these multi_hit recoreds
      multi_hit_genes = unique(multi_hits[, unlist(all_genes)])
      if (cesa@maf[, .N] > 0) {
        non_intergenic = cesa@maf[variant_type == "snv" & ! variant_id %in% cesa@mutations$snv[intergenic == T, snv_id]]
        maf_gene_counts = sort(table(unlist(non_intergenic$genes))[multi_hit_genes], decreasing = T)
      } else {
        maf_gene_counts = table(NA)
      }
      chosen_aac = character(num_to_check)
      
      # poorly optimized; may need to improve if large data sets are analyzed with lots of
      # overlapping AACs
      original_multi = copy(multi_hits)
      multi_hits = multi_hits[, .SD, .SDcols = c("variant_id", "aa_ref", "aa_alt", "essential_splice", "gene", "all_aac", maf_freq_cols)]
      all_aac_list = multi_hits$all_aac
      for (i in 1:num_to_check) {
        # if already picked an AAC out of this group, skip
        curr_aac = all_aac_list[[i]]
        if (any(curr_aac %in% chosen_aac)) {
          next
        }
        candidates = multi_hits[curr_aac, nomatch = NULL]
        multi_hits = multi_hits[! candidates$variant_id]

        # if other candidates are not in table for some reason, use the original entry
        if (candidates[, .N] == 1) {
          chosen_aac[i] = candidates$variant_id
          next
        }
        
        # take highest MAF frequency first (will usually be the same for all)
        if (cesa@maf[, .N] > 0) {
          freq_col = maf_freq_cols[length(maf_freq_cols)] # last freq col gives total frequency
          highest_freq = max(candidates[[freq_col]])
          candidates = candidates[candidates[[freq_col]] == highest_freq]
        }
        
        # then, prioritize splice-site variants
        if (any(candidates$essential_splice)) {
          candidates = candidates[essential_splice == T]
        }
        
        # next, choose non-silent over silent
        if (any(candidates$aa_ref != candidates$aa_alt)) {
          candidates = candidates[aa_ref != aa_alt]
        }
        # next, choose the variant(s) with the most mutations in the gene
        # If MAF is empty or has no gene records, skip
        if (length(maf_gene_counts) > 0) {
          candidates[, gene_var_count := as.numeric(maf_gene_counts[gene])]
          # Will be NA if variants are not in MAF and there are also no same-gene records in MAF
          if(! all(is.na(candidates$gene_var_count))) {
            max_count = max(candidates$gene_var_count, na.rm = T)
            candidates = candidates[gene_var_count == max_count]
          }
        }

        # arbitrarily return first alphabetically of remaining
        chosen_aac[i] = sort(candidates$variant_id)[1]
      }
      multi_hits = original_multi
      # remove secondary (non-chosen) AACs, but save all SNV IDs and re-select those passing filters
      chosen_aac = chosen_aac[chosen_aac != ""]
      not_chosen_aac = multi_hits[! chosen_aac, variant_id]
      all_const_snv = unlist(multi_hits$constituent_snvs)
      remaining_const_snv = multi_hits[chosen_aac, unlist(constituent_snvs)]
      combined = combined[! not_chosen_aac, on = "variant_id"]
      
      # Need to get annotations for SNVs that passed user's filters (given in subvariant_snv_ids),
      # and that are now no longer constituent SNVs
      snv_to_reselect = intersect(snvs_to_save, setdiff(all_const_snv, remaining_const_snv))
      
      if (length(snv_to_reselect) > 0) {
        reselected = select_variants(cesa, variant_passlist = snv_to_reselect)
        combined = rbind(combined, reselected[snv_to_reselect, on = "variant_id"])
      }
    }
  }
  
  # confirm that output table will be eligible for functions that require non-overlapping variants (e.g., ces_variant)
  # presume overlapping when remove_secondary_aac == FALSE
  nonoverlapping = FALSE 
  if (remove_secondary_aac) {
    if (length(intersect(combined[variant_type == "snv", variant_id], unlist(combined[variant_type == "aac", constituent_snvs]))) == 0) {
      nonoverlapping = TRUE
    } else {
      pretty_message(paste0("FYI, your output has overlapping variants (as in, it lists both an amino-acid-change variant and a constituent SNV as separate records). ",
                     "That means the output can't be fed into functions like ces_variant() that assume non-overlapping variants."))
    }
  }
  

  # order output in chr/pos order
  setkey(combined, "chr")
  combined = combined[BSgenome::seqnames(bsg), nomatch = NULL][, .SD[order(start)], by = "chr"]
  setcolorder(combined, c("variant_name", "variant_type", "chr", "start", "end", "variant_id", "ref", "alt", "gene", 
                          "strand", "aachange", "essential_splice", "intergenic", "trinuc_mut", "aa_ref", "aa_pos", "aa_alt", "coding_seq", 
                          "center_nt_pos", "pid", "constituent_snvs", "multi_anno_site", "all_aac", "all_genes",
                          "covered_in", unlist(zipup(all_cov_cols, maf_freq_cols))))
  
  setattr(combined, "cesa_id", cesa@advanced$uid)
  setattr(combined, "nonoverlapping", nonoverlapping)
  
  # May want to keep track of what covered_regions are in the CESAnalysis, since it's
  # possible to add more at any time. However, it's not possible to change covered_regions
  # of already-loaded samples, so plugging an out-of-data variant table into ces_variant
  # shouldn't ever affect output.
  #setattr(combined, "cesa_cov_regions", unlist(sapply(cesa@coverage, names), use.names = F))
  return(combined[]) # brackets force the output to print when unassigned (should automatically, but this is a known data.table issue)
}


#' Add variant annotations
#' 
#' Use this function to add variant annotations to your CESAnalysis by specifying variants
#' to add in one of five ways: a data.table containing genomic coordinates (output from
#' select_variants(), typically), a GRanges object, a BED file, another CESAnalysis, or
#' SNV IDs.
#' 
#' All methods of adding variants work by identifying which SNVs to add and then using the
#' target_cesa's associated reference data to identify overlapping amino-acid-change
#' mutations, which are then added as well. (You can't add just SNVs or just AACs.) Note
#' that if you try to add far more distinct variants than appear in a typical cohort (as
#' in, millions), annotation will take a while and the annotation tables in the
#' CESAnalysis may take up significant memory. Please contact us if you have issues.
#' 
#' @param target_cesa CESAnalysis to receive variant annotations
#' @param variant_table A data.table with chr/start/end positions (1-based closed
#'   coordinates, like MAF format). All possible SNVs overlapping the table's genomic
#'   coordinates (within \code{padding} bases) will be added. The tables returned by
#'   select_variants() and (CESAnalysis)$variants work, and get special handling of
#'   amino-acid-change SNVs: only the precise positions in start, end, and center_nt_pos
#'   are used. (This avoids adding all variants between start/end, which on
#'   splice-site-spanning variants can be many thousands.)
#' @param bed A path to a BED file. All possible SNVs overlapping BED intervals (within
#'   \code{padding} bases) will be added.
#' @param gr A GRanges object. All possible SNVs overlapping the ranges (within \code{padding}
#'   bases) will be added.
#' @param snv_id Character vector of CES-style SNV IDs to add.
#' @param source_cesa Another CESAnalysis from which to copy snv_ids. SNVs will be
#'   re-annotated using the target_cesa's associated reference data.
#' @param padding How many bases (default 0) to expand start and end of each gr range
#' @export
add_variants = function(target_cesa = NULL, variant_table = NULL, snv_id = NULL, bed = NULL, 
                        gr = NULL, source_cesa = NULL, padding = 0) {
  if(! is(target_cesa, "CESAnalysis")) {
    stop("target_cesa should be a CESAnalysis", call. = F)
  }
  target_cesa = update_cesa_history(target_cesa, match.call())
  if(! identical(target_cesa@advanced$annotated, T)) {
    # It's okay if there's no data in the CESAnalysis
    if(target_cesa@maf[, .N] == 0) {
      target_cesa@advanced$annotated = T
    } else {
      stop("target_cesa should be annotated", call. = F)
    }
  }

  if (! is(padding, "numeric") || length(padding) != 1 || trunc(padding) != padding || padding < 0) {
    stop("padding should be 1-length non-negative integer")
  }
  
  # just one possible method should be chosen
  if (sum(sapply(list(variant_table, gr, bed, snv_id, source_cesa), is.null)) != 4){
    stop("Exactly one method of adding variants must be chosen.")
  }


  if (! is.null(source_cesa)) {
    if(! is(source_cesa, "CESAnalysis")) {
      stop("source_cesa should be a CESAnalysis", call. = F)
    }
    source_snv_table = source_cesa@mutations$snv
    if (is.null(source_snv_table)) {
      stop("source_cesa has no SNV annotations", call. = F)
    }
    
    if (! identical(target_cesa@ref_key, source_cesa@ref_key)) {
      stop("The pair of CESAnalysis objects appear to use different reference data sets.")
    } else {
      if (! target_cesa@ref_key %in% names(.official_refsets)) {
        if(! identical(target_cesa@ref_data_dir, source_cesa@ref_data_dir)) {
          msg = paste0("Custom refset data directories may differ (", target_cesa@ref_data_dir, 
                 ", ", source_cesa@ref_data_dir, "). If the refsets are not equivalent, annotations will be corrupted.")
          warning(pretty_message(msg, emit = F))
        }
      }
      if(length(target_cesa@mutations) == 0) {
        target_cesa@mutations$snv = copy(source_cesa@mutations$snv)
        target_cesa@mutations$amino_acid_change = copy(source_cesa@mutations$amino_acid_change)
      } else {
        new_snvs = source_cesa@mutations$snv[! target_cesa@mutations$snv$snv_id, on = "snv_id"]
        if (new_snvs[, .N] == 0) {
          stop("There are no new variants to copy over.")
        }
        # covered_in may vary, but doesn't matter because it will be regenerated from scratch
        target_cesa@mutations$snv = rbind(target_cesa@mutations$snv, new_snvs)
        new_aacs = source_cesa@mutations$snv[! target_cesa@mutations$amino_acid_change$aac_id, on = "aac_id"]
        target_cesa@mutations$amino_acid_change = rbind(target_cesa@mutations$amino_acid_change, new_aacs)
      }
      return(update_covered_in(target_cesa))
    }
  }
  
  # Handle gr, bed, variant_table: all get converted to a validated gr before creation of SNV IDs
  # We've already ensured that only one of these can be non-null
  if (! all(sapply(list(variant_table, gr, bed), is.null))) {
    bsg = get_cesa_bsg(target_cesa)
    if (! is.null(bed)) {
      if(is.character(bed)) {
        if (length(bed) != 1) {
          stop("bed should be a path (one-length character) to a BED file.")
        }
        if (! file.exists(bed)) {
          stop("BED file not found (check path?)")
        }
        gr = rtracklayer::import.bed(bed)
      } else {
        stop("bed should be a path (one-length character) to a BED file.")
      }
    }
    if (! is.null(gr)) {
      if (! is(gr, "GRanges")) {
        stop("gr should be a GRanges object")
      }
    }
    if (! is.null(variant_table)) {
      gr = get_gr_from_table(variant_table)
    }
    
    # Validate GRanges and add padding if specified
    gr = clean_granges_for_cesa(cesa = target_cesa, gr = gr, padding = padding)
    if(length(gr) == 0) {
      stop("No variants present in the input.")
    }
    # convert to GPos and put in MAF-like table
    gpos = GenomicRanges::GPos(gr)
    ref = as.character(BSgenome::getSeq(bsg, gpos))
    snv_table = data.table(chr = as.character(GenomicRanges::seqnames(gpos)), pos = GenomicRanges::pos(gpos),
                           ref = ref)
    nt = c("A", "C", "G", "T")
    if (any(! ref %in% nt)) {
      snv_table = snv_table[ref %in% nt]
      message("Note: some variants in input were dropped because of ambiguous reference sequence (N's)")
    }
    snv_table = snv_table[rep(1:.N, each = 4)]
    snv_table[, alt := rep.int(c("A", "C", "G", "T"), .N/4)]
    snv_table = snv_table[ref != alt]
    snvs_to_annotate = snv_table[, paste0(chr, ':', pos, '_', ref, '>', alt)]
  }
  
  # If supplied SNV IDs (rather than source_cesa, gr, bed, variant_table), validate them
  if(! is.null(snv_id)) {
    if(! is(snv_id, "character") | length(snv_id) == 0) {
      stop("Expected snv_id to be character vector of snv_ids (e.g., 1:100_A>G", call. = F)
    }
    # will stop with errror if any IDs fail validation
    validate_snv_ids(snv_id, get_cesa_bsg(target_cesa))
    snvs_to_annotate = snv_id
  }
  
  num_variants = length(snvs_to_annotate)
  if(num_variants == 0) {
    stop("No SNVs to add (check your input).")
  }
  snvs_to_annotate = setdiff(snvs_to_annotate, target_cesa@mutations$snv$snv_id)
  num_to_add = length(snvs_to_annotate)
  
  if(num_to_add == 0) {
    stop("Tried to add ", num_variants , " variants, but all of them are already annotated in the CESAnalysis.")
  }
  
  num_identified_str = format(num_variants, big.mark = ",")
  num_to_add_str = format(num_to_add, big.mark = ",")
  
  if (num_to_add == num_variants) {
    pretty_message(paste0("Annotating ", num_to_add_str, " variants..."))
  } else {
    pretty_message(paste0("Received ", num_identified_str, " total variants ", 
                          " and annotating the ", num_to_add_str, " variants that are new..."))
  }
  
  if(num_to_add > 1e6) {
    warning("You're adding a lot of variants! Let us know if you have any issues.", immediate. = T, call. = F)
  }
  
  # convert snv_id into MAF-like table, with an NA UPI
  cesa = target_cesa
  snv_id = snvs_to_annotate
  maf = as.data.table(tstrsplit(snv_id, split = '[:_>]'))
  colnames(maf) = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
  maf$Unique_Patient_Identifier = NA_character_
  setcolorder(maf, "Unique_Patient_Identifier")
  maf$variant_type = "snv"
  maf[, Start_Position := as.numeric(Start_Position)]
  
  # add the new variants to cesa@maf (will remove after annotation)
  cesa@maf = rbind(cesa@maf, maf, fill = T)
  prev_recording_status = cesa@advanced$recording
  cesa@advanced$recording = F
  cesa = annotate_variants(cesa)
  cesa@advanced$recording = prev_recording_status
  cesa@maf = cesa@maf[! is.na(Unique_Patient_Identifier)]
  return(cesa)
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
  
  if (! all(grepl('^[1-9][0-9]*$', dt$start))) {
    stop("Some SNV IDs have illegal positions (watch out for mistakes like \"2:1e+06_A>C\").")
  }
  
  
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

#' Get GRanges from chr/start/end table
#' 
#' Mainly built for select_variants() output, and uses the
#' center_nt_pos on AACs (rather than all from start-end). Assumes
#' MAF-like coordinates (1-based, closed).
#' 
#' @param variant_table data.table
#' @keywords internal
get_gr_from_table = function(variant_table) {
  if (! is(variant_table, "data.table")) {
    stop("variant_table should be a data.table.")
  }
  
  if (! all(c("chr", "start", "end") %in% names(variant_table))) {
    stop("variant table should have chr/start/end columns.")
  }
  
  # Use just the start, end, and center_nt_pos positions (not all positions from start to end)
  if (all(c("variant_type", "center_nt_pos") %in% names(variant_table))) {
    non_aac = variant_table[variant_type != 'aac', .(chr, start, end)]
    aac_table = variant_table[variant_type == 'aac']
    aac_table = rbindlist(list(aac_table[, .(chr, start, end = start)], 
                               aac_table[, .(chr, start = end, end)],
                               aac_table[, .(chr, start = center_nt_pos, end = center_nt_pos)]))
    variant_table = rbind(non_aac, aac_table)
  }
  gr = GenomicRanges::makeGRangesFromDataFrame(variant_table, ignore.strand = T)
  return(gr)
}






