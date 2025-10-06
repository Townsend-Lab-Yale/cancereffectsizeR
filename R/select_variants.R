#' Select and filter variants
#'
#' This function helps you find and view variant data from your CESAnalysis's MAF data and
#' mutation annotation tables. By default, substitutions that occur within coding regions are represented at
#' the codon level as amino acid changes (including synonymous substitutions, so amino acid "change"
#' is a bit of a misnomer), while substitutions outside of coding regions are represented as SBS. You can apply a series of filters to restrict output to
#' certain genes or genomic regions, require a minimum variant frequency in MAF data, and/or
#' specify exact \code{variant_ids} to return.
#' 
#' Only variants that are annotated in the CESAnalysis can be returned. To view more annotations,
#' such as for variants absent from the MAF data, you must first call add_variants() to add them to
#' the CESAnalysis.
#' 
#' Note that while intergenic SBS have their nearest genes annotated in the SBS
#' tables, these variants will not be captured by the genes filter of this function.
#' 
#' Definitions of some less self-explanatory columns:
#' \itemize{
#'   \item variant_name: In a coding variant, gene and protein change on the (MANE) canonical
#'   transcripts, such as "BRAF V600E". For coding changes reported on other transcripts, the
#'   protein ID is included: "POLH W415C (ENSP00000361300.1)". With older reference data sets
#'   (ces.refset.hg19, versions of ces.refset.hg38 < 1.3, and any custom reference data set that
#'   doesn't have complete information on canonical transcripts), the variant name is a 
#'   shortening of the variant_id.
#'   \item start/end: lowest/highest genomic positions overlapping variant
#'   \item variant_id: unique IDs for variants given the associated genome assembly version and the transcript data
#'   \item ref/alt: genomic reference and alternate alleles (for genomic
#'   variants; NA for AACs) 
#'   \item gene: the affected gene in AACs; for SBS, the
#'   overlapping gene (or an arbitrary gene if more than one overlaps), or the nearest gene
#'   for intergenic SBS
#'   \item strand: for AACs, 1 if the reference sequence strand is the coding strand; -1 otherwise
#'   \item essential_splice: Variant is 1,2 bp upstream or
#'   1,2,5 bp downstream of an annotated splice position (edge case: if an SBS has
#'   multiple gene/transcript annotations, this doesn't say which one it's essential for)
#'   \item intergenic: variant does not overlap any coding regions in the reference data
#'   \item trinuc_mut: for SBS, the reference trinucleotide context, in deconstructSigs notation
#'   \item coding_seq: coding strand nucleotides in order of transcription
#'   \item center_nt_pos: regardless of strand, start/end give positions of two out of three AAC nucleotides; this
#'                       gives the position of the center nucleotide (maybe useful if the AAC spans a splice site)
#'   \item maf_prevalence: number of occurrences of the variant in MAF data
#'   \item samples_covering: number of MAF samples with sequencing coverage at the variant site
#' }
#' 
#' @param cesa CESAnalysis with MAF data loaded and annotated (e.g., with \code{load_maf()})
#' @param genes Filter variants to specified genes.
#' @param min_freq Filter out variants with MAF frequency below threshold (default 0).
#'   Note that variants that are not in the annotation tables will never be returned. Use
#'   \code{add_variants()} to include variants absent from MAF data in your CESAnalysis.
#' @param variant_ids Filter out any variants besides those specified. (Variants specified here can
#'   still get removed by other filters.)
#' @param type Filter out variants not of the specified types. Supported types are sbs, aac, dbs, dbs_aac.
#' @param gr Filter out any variants not within input GRanges +/- \code{padding} bases.
#' @param variant_position_table Filter out any variants that don't intersect the
#'   positions given in chr/start/end of this table (1-based closed coordinates).
#'   Typically, the table comes from a previous \code{select_variants} call and can be expanded
#'   with \code{padding}. (Gritty detail: Amino acid substitutions get special handling. Only the
#'   precise positions in start, end, and center_nt_pos are used. Otherwise, coding changes that
#'   span splice sites would capture all the intronic sites.)
#' @param padding add +/- this many bases to every range specified in \code{gr} or
#'   \code{variant_position_table} (stopping at chromosome ends, naturally).
#' @return A data.table with info on selected variants (see details)
#' @export
select_variants = function(cesa, genes = NULL, min_freq = 0, variant_ids = NULL, type = NULL, gr = NULL, variant_position_table = NULL, 
                            padding = 0) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis object")
  }
  bsg = get_cesa_bsg(cesa)
  
  # Soon we'll check indels, too
  if(cesa@mutations$sbs[, .N] + cesa@mutations$dbs[, .N] == 0) {
    stop("There are no variants in the CESAnalysis.")
  }
  
  if(! is.numeric(min_freq) || length(min_freq) > 1 || min_freq - as.integer(min_freq) != 0) {
    stop("min_freq should be 1-length integer at least zero")
  }
  if(! is.null(genes) && ! is(genes, "character")) {
    stop("genes should be character vector of genes to include or NULL (for no gene filtering)")
  }
  
  if(! is.null(variant_ids) && ! is(variant_ids, "character")) {
    stop("variant_ids should be character vector of variant IDs to include, or left NULL.")
  }
  
  if(! is.null(type)) {
    if(! is.character(type)) {
      stop('type must be a character vector of variant types.')
    }
    type = unique(tolower(type))
    if(length(setdiff(type, c('sbs', 'aac', 'dbs', 'dbs_aac'))) > 0) {
      msg = "Unsupported variant types specified."
      stop(pretty_message(msg, emit = F))
    }
  } else {
    type = c('sbs', 'aac', 'dbs', 'dbs_aac')
  }
  
  # Start with all variants, then apply filters
  
  selected_aac_ids = character()
  if('aac' %in% type) {
    selected_aac_ids = cesa@mutations$amino_acid_change$aac_id
  }
  selected_sbs_ids = character()
  if('sbs' %in% type) {
    selected_sbs_ids = cesa@mutations$sbs$sbs_id
  }
  selected_dbs_ids = character()
  if('dbs' %in% type) {
    selected_dbs_ids = cesa@mutations$dbs$dbs_id
  }
  selected_dbs_aac_ids = character()
  if('dbs_aac' %in% type) {
    selected_dbs_aac_ids = cesa@mutations$dbs_codon_change$dbs_aac_id
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
    sbs_gr = GenomicRanges::makeGRangesFromDataFrame(cesa@mutations$sbs, seqnames.field = "chr", start.field = "pos", 
                                                           end.field = "pos", seqinfo = GenomeInfoDb::seqinfo(final_gr))
    captured_sbs = cesa@mutations$sbs[IRanges::overlapsAny(query = sbs_gr, subject = final_gr, type = "within")]
    
    all_dbs = cesa@mutations$dbs[, .(chr, start = pos, end = pos + 1)]
    dbs_gr = GenomicRanges::makeGRangesFromDataFrame(all_dbs, seqinfo = GenomeInfoDb::seqinfo(final_gr))
    captured_dbs = cesa@mutations$dbs[IRanges::overlapsAny(query = dbs_gr, subject = final_gr, type = "within")]
    
    if (captured_sbs[, .N] + captured_dbs[, .N] == 0) {
      return(data.table())
    }
    selected_sbs_ids = intersect(selected_sbs_ids, captured_sbs$sbs_id)
    aac_passing_gr = cesa@mutations$aac_sbs_key[captured_sbs$sbs_id, unique(aac_id), on = 'sbs_id', nomatch = NULL]
    selected_aac_ids = intersect(selected_aac_ids, aac_passing_gr)
    
    selected_dbs_ids = intersect(selected_dbs_ids, captured_dbs$dbs_id)
    dbs_aac_passing_gr = cesa@mutations$aac_dbs_key[captured_dbs$dbs_id, unique(dbs_aac_id), on = 'dbs_id', nomatch = NULL]
    selected_dbs_aac_ids = intersect(selected_dbs_aac_ids, dbs_aac_passing_gr)
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
    # Note we're not returning intergenic sbs that just have one of the chosen genes as their nearest gene
    genes_by_sbs = cesa@mutations$sbs[intergenic == FALSE, .(gene = unlist(genes)), by = "sbs_id"]
    sbs_in_genes = genes_by_sbs[gene %in% genes, unique(sbs_id)]
    selected_sbs_ids = intersect(selected_sbs_ids, sbs_in_genes)
    
    dbs_aac_in_genes = cesa@mutations$dbs_codon_change[gene %in% genes, dbs_aac_id]
    selected_dbs_aac_ids = intersect(selected_dbs_aac_ids, dbs_aac_in_genes)
    
    # KNOWN ISSUE: No nearest gene information for DBS
    selected_dbs_ids = intersect(selected_dbs_ids, cesa@mutations$aac_dbs_key[dbs_aac_id %in% dbs_aac_in_genes, dbs_id])
  }
  
  # Filter by variant_id (in other words, discard all other variants)
  if (length(variant_ids) > 0) {
    variant_ids = unique(variant_ids)
    matching_sbs_ids = cesa@mutations$sbs[variant_ids, sbs_id, nomatch = NULL]
    matching_aac_ids = cesa@mutations$amino_acid_change[variant_ids, aac_id, nomatch = NULL]
    matching_dbs_ids = cesa@mutations$dbs[variant_ids, dbs_id, on = 'dbs_id', nomatch = NULL]
    matching_dbs_aac_ids = cesa@mutations$dbs_codon_change[variant_ids, dbs_aac_id, on = 'dbs_aac_id', nomatch = NULL]
    
    # KNOWN ISSUE: spaces should be accepted in SBS/DBS ID
    missing_ids = setdiff(variant_ids, c(matching_sbs_ids, matching_aac_ids, matching_dbs_ids, matching_dbs_aac_ids))
    
    # if any IDs are missing, try to interpret them as "short" AAC names (i.e., without protein ID)
    if (length(missing_ids) > 0) {
      aac_matches = cesa@mutations$amino_acid_change[missing_ids, on = "variant_name"]
      
      missing_ids = aac_matches[is.na(aac_id), variant_name]
      
      dbs_aac_matches = cesa@mutations$dbs_codon_change[missing_ids, on = 'variant_name']
      missing_ids = dbs_aac_matches[is.na(dbs_aac_id), variant_name]
      
      if (length(missing_ids) > 0) {
        num_missing = length(missing_ids)
        missing_ids = paste(missing_ids, collapse = ", ")
        stop(num_missing, " variants given in variant_ids couldn't be found. Either they're bad IDs or there are no annotations for them.\n",
             missing_ids)
      }
      matching_aac_ids = c(matching_aac_ids, aac_matches$aac_id)
      matching_dbs_aac_ids = c(matching_dbs_aac_ids, dbs_aac_matches$dbs_aac_id)
      if (any(duplicated(aac_matches$variant_name) || any(duplicated(dbs_aac_matches$variant_name)))) {
        msg = paste0("Amino acid changes with no protein ID (styled like \"KRAS G12C\") were recognized and matched ",
                     "with specific variant IDs. However, some of your variant names matched more than ",
                     "one aac_id (i.e., the same amino acid change is possible on multiple protein isoforms, and both are present ",
                     "in this analysis's variant annotations.) When the underlying SBS are the same, only one AAC will be returned ",
                     "Otherwise, all matching variants will be returned. ",
                     "(There may be additional matching variants with MAF frequency = 0 that are not ",
                     "annotated in this analysis; these will not be returned.)")
        pretty_message(msg, black = F)
      }
    }
    selected_sbs_ids = intersect(selected_sbs_ids, matching_sbs_ids)
    selected_aac_ids = intersect(selected_aac_ids, matching_aac_ids)
    selected_dbs_ids = intersect(selected_dbs_ids, matching_dbs_ids)
    selected_dbs_aac_ids = intersect(selected_dbs_aac_ids, matching_dbs_aac_ids)
  }
  
  
  # Constituent SBS of AACs that get tossed in de-overlapping need to be
  # saved if they pass filters in their own right and don't overlap the AACs kept. Keep
  # track here of SBS that need to appear in final output. Note if any of these don't
  # pass frequency filter, they shouldn't be saved, and they get removed below.
  sbs_to_recover = selected_sbs_ids
  dbs_to_recover = selected_dbs_ids
  
  # Remove variants contained in other variants, unless variants were explicitly specified by variant_id
  aac_sbs_key = cesa@mutations$aac_sbs_key[selected_aac_ids, on = 'aac_id']
  aac_dbs_key = cesa@mutations$aac_dbs_key[selected_dbs_aac_ids, on = 'dbs_aac_id']
  
  selected_sbs_ids = setdiff(selected_sbs_ids, aac_sbs_key$sbs_id)
  selected_dbs_ids = setdiff(selected_dbs_ids, aac_dbs_key$dbs_id)
  
  selected_sbs = setDT(cesa@mutations$sbs[selected_sbs_ids, on = 'sbs_id'])
  selected_aac = setDT(cesa@mutations$amino_acid_change[selected_aac_ids, on = 'aac_id'])
  selected_dbs = setDT(cesa@mutations$dbs[selected_dbs_ids, on = 'dbs_id'])
  selected_dbs_aac = setDT(cesa@mutations$dbs_codon_change[selected_dbs_aac_ids, on = 'dbs_aac_id'])
  
  # Get variant counts and coverage
  sbs_from_aac = cesa@mutations$aac_sbs_key[selected_aac$aac_id, .(aac_id, sbs_id), on = 'aac_id']
  dbs_from_aac = cesa@mutations$aac_dbs_key[selected_dbs_aac$dbs_aac_id, .(dbs_aac_id, dbs_id), on = 'dbs_aac_id']
  
  if(sbs_from_aac[, .N] == 0 && length(selected_sbs_ids) == 0  && dbs_from_aac[, .N] == 0 && length(selected_dbs_ids) == 0) {
    message("No variants passed selection criteria!")
    return(NULL)
  }
  if (cesa@maf[, .N] > 0) {
    counts_and_cov = .variant_counts(cesa, samples = cesa@samples[, .(patient_id, covered_regions)],
                                     sbs_from_aac = sbs_from_aac,
                                     noncoding_sbs_id = selected_sbs_ids,
                                     dbs_from_aac = dbs_from_aac,
                                     noncoding_dbs_id = selected_dbs_ids)
    setnames(counts_and_cov, c("N", "num_cov"), c("maf_prevalence", "samples_covering"))
    
  } else {
    counts_and_cov = data.table(variant_id = c(selected_aac$aac_id, selected_sbs_ids),
                                variant_type = c(rep('aac', selected_aac[, .N]),
                                                 rep('sbs', length(selected_sbs_ids)),
                                                 rep('dbs_aac', selected_dbs_aac[, .N]),
                                                 rep('dbs', length(selected_dbs_ids))),
                                maf_prevalence = 0, samples_covering = 0)
  }
  
  selected_sbs[counts_and_cov, c("maf_prevalence", "samples_covering") := list(maf_prevalence, samples_covering), on = c(sbs_id = 'variant_id')]
  selected_sbs = selected_sbs[maf_prevalence >= min_freq]
  selected_aac[counts_and_cov, c("maf_prevalence", "samples_covering") := list(maf_prevalence, samples_covering), on = c(aac_id = 'variant_id')]
  selected_aac = selected_aac[maf_prevalence >= min_freq]
  selected_dbs[counts_and_cov, let(maf_prevalence = maf_prevalence, samples_covering = samples_covering),
               on = c(dbs_id = 'variant_id')]
  selected_dbs_aac[counts_and_cov, let(maf_prevalence = maf_prevalence, samples_covering = samples_covering),
               on = c(dbs_aac_id = 'variant_id')]
  
  
  combined = data.table()
  if(selected_sbs[, .N] > 0) {
    selected_sbs[, variant_type := "sbs"]
    selected_sbs[, strand := NA_integer_] # because AAC table is +1/-1
    selected_sbs[, c("start", "end") := .(pos, pos)]
    selected_sbs[, pos := NULL]
    # for convenience, take 1 gene per sbs for the output gene column
    selected_sbs[, gene := sapply(genes, function(x) x[1])] 
    selected_sbs[, c('genes', 'nearest_pid')] = NULL # not keeping these list-type annotations
    setnames(selected_sbs, "sbs_id", "variant_id")
    combined = rbind(combined, selected_sbs, fill = TRUE, use.names = TRUE)
  }
  
  if(selected_dbs[, .N] > 0) {
    # KNOWN ISSUE: no gene
    selected_dbs[, let(variant_type = 'dbs', variant_name = sub('_', ' ', dbs_id),
                       strand = NA_integer_, start = pos, end = pos + 1, pos = NULL,
                       gene = NA_character_)]
    setnames(selected_dbs, 'dbs_id', 'variant_id')
    combined = rbind(combined, selected_dbs, fill = TRUE, use.names = TRUE)
  }
  
  if(selected_aac[, .N] > 0) {
    selected_aac[, variant_type := "aac"]
    selected_aac[, intergenic := FALSE]
    selected_aac[, start := pmin(nt1_pos, nt3_pos)]
    selected_aac[, end := pmax(nt1_pos, nt3_pos)]
    selected_aac[, center_nt_pos := nt2_pos]
    selected_aac[, c("nt1_pos", "nt2_pos", "nt3_pos") := NULL]
    setnames(selected_aac, "aac_id", "variant_id")
    combined = rbind(combined, selected_aac, fill = TRUE, use.names = T)
  }
  
  if(selected_dbs_aac[, .N] > 0) {
    # KNOWN ISSUE: start/end NA
    selected_dbs_aac[, let(variant_type = 'dbs_aac', intergenic = FALSE,
                           start = NA, end = NA)]
    setnames(selected_dbs_aac, 'dbs_aac_id', 'variant_id')
    combined = rbind(combined, selected_dbs_aac, fill = TRUE, use.names = T)
  }
  
  # Combine SBS and AAC tables
  if(combined[, .N] == 0) {
    return(NULL)
  }
  
  # Handle overlapping  mutations using tiebreakers explained below.
  multi_anno_sbs = aac_sbs_key[, .N, by = 'sbs_id'][N > 1, sbs_id]
  multi_hit_aac = aac_sbs_key[sbs_id %in% multi_anno_sbs, unique(aac_id)]
  multi_anno_dbs = aac_dbs_key[, .N, by = 'dbs_id'][N > 1, dbs_id]
  multi_hit_dbs_aac = aac_dbs_key[dbs_id %in% multi_anno_dbs, unique(dbs_aac_id)]
  multi_hits = combined[variant_id %in% c(multi_hit_aac, multi_hit_dbs_aac)]
  aa_key = rbind(aac_sbs_key[aac_id %in% multi_hits$variant_id, .(nt_id = sbs_id, aa_id = aac_id)],
                 aac_dbs_key[dbs_aac_id %in% multi_hits$variant_id, .(nt_id = dbs_id, aa_id = dbs_aac_id)])
  multi_hits = merge.data.table(multi_hits, aa_key, by.x = 'variant_id', by.y = 'aa_id', all = TRUE)
  num_to_check = multi_hits[, .N]

  if (num_to_check > 0) {
    multi_hits[, is_premature := aa_alt == "STOP" & aa_ref != "STOP"]
        
    # Apply prioritization
    if(check_for_ref_data(cesa, 'transcripts')) {
        transcripts = unique(get_ref_data(cesa, 'transcripts')[, .(is_mane, is_mane_plus, pid = protein_id)])
        multi_hits[transcripts, c('is_mane', 'is_mane_plus') := .(is_mane, is_mane_plus), on = 'pid']
        multi_hits = multi_hits[order(-essential_splice, -is_premature, aa_ref == aa_alt, -is_mane, -is_mane_plus, -maf_prevalence, variant_id)]
        multi_hits[, c('is_mane', 'is_mane_plus') := NULL]
    } else {
        multi_hits = multi_hits[order(-essential_splice, -is_premature, aa_ref == aa_alt, -maf_prevalence, variant_id)]
    }
    multi_hits[, is_premature := NULL]
    multi_hits[, nt_already_used := duplicated(nt_id), by = 'variant_type'] # same nt may be in DBS AACs and SBS AACs
    multi_hits[, to_use := ! any(nt_already_used), by = 'variant_id']
    
    chosen_aa_variants = multi_hits[to_use == TRUE, variant_id]
    
    # remove secondary (non-chosen) AACs/DBS AACs, but save all SBS/DBS IDs and re-select those passing filters
    not_chosen_aa_variants = setdiff(multi_hits$variant_id, chosen_aa_variants)
    
    combined = combined[! not_chosen_aa_variants, on = 'variant_id']
    
    not_used = setdiff(multi_hits[to_use == F, nt_id], multi_hits[to_use == T, nt_id])
    to_reselect = intersect(not_used, c(sbs_to_recover, dbs_to_recover))
    
    # We need to recover sbs that were in AACs and that passed user's filters,
    # but that are now no longer constituent sbs after de-overlapping.
    if (length(to_reselect) > 0) {
      reselected = select_variants(cesa, variant_ids = to_reselect)
      reselected = reselected[, .SD, .SDcols = intersect(names(combined), names(reselected))]
      combined = rbind(combined, reselected, fill = TRUE)
    }
  }
  # order output in chr/pos order
  setkey(combined, "chr")
  combined = setDT(combined[order(start)][BSgenome::seqnames(bsg), nomatch = NULL, on = "chr"])
  
  # Add fields that may be missing depending on what variant types were selected.
  expected_fields = c("variant_name", "variant_type", "chr", "start", "end", "variant_id", "ref", "alt", "gene", 
                      "strand", "aachange", "essential_splice", "intergenic", "trinuc_mut", "cosmic_dbs_class", 
                      "aa_ref", "aa_pos", "aa_alt", "coding_seq", 
                      "center_nt_pos", "pid", "maf_prevalence", "samples_covering")
  missing_fields = setdiff(expected_fields, names(combined))
  if(length(missing_fields) > 0) {
    combined[, (missing_fields) := NA]
  }
  combined = combined[, .SD, .SDcols = expected_fields] # sometimes will have extra fields from reselecting
  
  if(check_for_ref_data(cesa, "transcripts")) {
    transcripts = unique(get_ref_data(cesa, "transcripts")[, .(is_mane, transcript_tags, pid = protein_id)])
    combined[transcripts, let(is_MANE_transcript = is_mane,
                              transcript_tags = transcript_tags),
                              on = 'pid']
    combined[variant_type %in% c('sbs', 'dbs'), is_MANE_transcript := NA]
  }
  setcolorder(combined, expected_fields)

  setattr(combined, "cesa_id", cesa@advanced$uid)
  setkey(combined, 'variant_id', physical = F)
  return(combined[]) # brackets force the output to print when unassigned (should automatically, but this is a known data.table issue)
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






