#' Select and filter variants
#'
#' This function helps you find and view variant data from your CESAnalysis's MAF data and
#' mutation annotation tables. By default, almost all amino-acid-changing mutations and
#' noncoding SBS are returned. You can apply a series of filters to restrict output to
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
#'   \item multi_anno_site: T/F whether variant has multiple gene/transcript/AAC annotations
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
#' @param type Vector of variant types to include, such as "sbs" or "aac".
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
    if(length(setdiff(type, c('sbs', 'aac'))) > 0) {
      msg = "Currently, annotations are only available for SBS (single-base substitution) and AAC (amino-acid-change) variant types."
      stop(pretty_message(msg, emit = F))
    }
  } else {
    type = c('sbs', 'aac')
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
    mutations_gr = GenomicRanges::makeGRangesFromDataFrame(cesa@mutations$sbs, seqnames.field = "chr", start.field = "pos", 
                                                           end.field = "pos", seqinfo = GenomeInfoDb::seqinfo(final_gr))
    captured = cesa@mutations$sbs[IRanges::overlapsAny(query = mutations_gr, subject = final_gr, type = "within")]
    if (captured[, .N] == 0) {
      return(data.table())
    }
    selected_sbs_ids = intersect(selected_sbs_ids, captured$sbs_id)
    aac_passing_gr = cesa@mutations$aac_sbs_key[captured$sbs_id, unique(aac_id), on = 'sbs_id', nomatch = NULL]
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
    
    # Note we're not returning intergenic sbs that just have one of the chosen genes as their nearest gene
    genes_by_sbs = cesa@mutations$sbs[intergenic == FALSE, .(gene = unlist(genes)), by = "sbs_id"]
    sbs_in_genes = genes_by_sbs[gene %in% genes, unique(sbs_id)]
    selected_sbs_ids = intersect(selected_sbs_ids, sbs_in_genes)
  }
  
  # Filter by variant_id (in other words, discard all other variants)
  if (length(variant_ids) > 0) {
    variant_ids = unique(variant_ids)
    matching_sbs_ids = cesa@mutations$sbs[variant_ids, sbs_id, nomatch = NULL]
    matching_aac_ids = cesa@mutations$amino_acid_change[variant_ids, aac_id, nomatch = NULL]
    missing_ids = setdiff(variant_ids, c(matching_sbs_ids, matching_aac_ids))
    
    # if any IDs are missing, try to interpret them as "short" AAC names (i.e., without protein ID)
    if (length(missing_ids) > 0) {
      tmp = cesa@mutations$amino_acid_change[, .(aac_id, variant_name)]
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
        msg = paste0("Shorthand amino-acid-change names (styled like \"KRAS_G12C\") were recognized and matched ",
                     "with cancereffectsizeR-style aac_ids. However, some of your variant names matched more than ",
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
  }
  
  
  
  # Constituent SBS of AACs that get tossed in de-overlapping need to be
  # saved if they pass filters in their own right and don't overlap the AACs kept. Keep
  # track here of SBS that need to appear in final output. Note if any of these don't
  # pass frequency filter, they shouldn't be saved, and they get removed below.
  sbs_to_recover = selected_sbs_ids
  
  # Remove variants contained in other variants, unless variants were explicitly specified by variant_id
  aac_sbs_key = cesa@mutations$aac_sbs_key[selected_aac_ids, on = 'aac_id']
  
  if(is.null(variant_ids)) {
    selected_sbs_ids = setdiff(selected_sbs_ids, aac_sbs_key$sbs_id)
  }
  
  selected_sbs = setDT(cesa@mutations$sbs[selected_sbs_ids, on = 'sbs_id'])
  selected_aac = setDT(cesa@mutations$amino_acid_change[selected_aac_ids, on = 'aac_id'])
  
  # Get variant counts and coverage
  sbs_from_aac = cesa@mutations$aac_sbs_key[selected_aac$aac_id, .(aac_id, sbs_id), on = 'aac_id']
  
  if(sbs_from_aac[, .N] == 0 && length(selected_sbs_ids) == 0) {
    message("No variants passed selection criteria!")
    return(NULL)
  }

  if (cesa@maf[, .N] > 0) {
    counts_and_cov = .variant_counts(cesa, samples = cesa@samples[, .(patient_id, covered_regions)],
                                     sbs_from_aac = aac_sbs_key[, .(aac_id, sbs_id)],
                                     noncoding_sbs_id = selected_sbs_ids)
    setnames(counts_and_cov, c("N", "num_cov"), c("maf_prevalence", "samples_covering"))
    
  } else {
    counts_and_cov = data.table(variant_id = c(selected_aac$aac_id, selected_sbs_ids),
                                variant_type = c(rep('aac', selected_aac[, .N]), rep('sbs', length(selected_sbs_ids))),
                                maf_prevalence = 0, samples_covering = 0)
  }
  
  selected_sbs[counts_and_cov, c("maf_prevalence", "samples_covering") := list(maf_prevalence, samples_covering), on = c(sbs_id = 'variant_id')]
  selected_sbs = selected_sbs[maf_prevalence >= min_freq]
  selected_aac[counts_and_cov, c("maf_prevalence", "samples_covering") := list(maf_prevalence, samples_covering), on = c(aac_id = 'variant_id')]
  selected_aac = selected_aac[maf_prevalence >= min_freq]
  
  # Annotate sbs table and prepare to merge with AACs
  selected_sbs[, variant_type := "sbs"]
  selected_sbs[, variant_name := sub('_', ' ', sbs_id)] # SBS IDs are already short and uniquely identifying
  selected_sbs[, strand := NA_integer_] # because AAC table is +1/-1
  selected_sbs[, c("start", "end") := .(pos, pos)]
  selected_sbs[, pos := NULL]
  # for convenience, take 1 gene per sbs for the output gene column
  if (selected_sbs[, .N] > 0) {
    selected_sbs[, gene := sapply(genes, function(x) x[1])] 
  } else {
    selected_sbs[, gene := character()]
  }
  selected_sbs[cesa@mutations$aac_sbs_key, multi_anno_site := multi_anno_site, on = 'sbs_id']
  selected_sbs[is.na(multi_anno_site), multi_anno_site := FALSE]
  selected_sbs[, c('genes', 'nearest_pid')] = NULL # not keeping these list-type annotations
  setnames(selected_sbs, "sbs_id", "variant_id")
  
  # AACs get a short variant name that might not be uniquely identifying if a gene has more than one CDS
  selected_aac[, variant_type := "aac"]
  selected_aac[, intergenic := FALSE]
  selected_aac[, start := pmin(nt1_pos, nt3_pos)]
  selected_aac[, end := pmax(nt1_pos, nt3_pos)]
  selected_aac[, center_nt_pos := nt2_pos]
  selected_aac[, c("nt1_pos", "nt2_pos", "nt3_pos") := NULL]
  
  # We call an AAC a multi-anno site if any of its SBS has multiple annotations.
  multi_anno_aac = cesa@mutations$aac_sbs_key[, any(multi_anno_site), by = 'aac_id']
  selected_aac[multi_anno_aac, multi_anno_site := V1, on = 'aac_id']
  setnames(selected_aac, "aac_id", "variant_id")
  
  # Combine SBS and AAC tables
  combined = rbindlist(list(selected_aac, selected_sbs), use.names = T, fill = T)
  if(combined[, .N] == 0) {
    message("No variants passed selection criteria!")
    return(NULL)
  }
  
  # Handle overlapping  mutations using tiebreakers explained below.
  # Exception: Variants specified by variant_id always get included.
  multi_hits = combined[variant_type == "aac" & multi_anno_site == TRUE]
  num_to_check = multi_hits[, .N]
  if (num_to_check > 0 && is.null(variant_ids)) {
    # for tie-breaking, count how many mutations are in each gene found in these multi_hit records
    # will need to produce these counts from scratch since some of the variants may not be in this select_variants() run
    aac_sbs_key = cesa@mutations$aac_sbs_key[multi_anno_site == TRUE] # only multi-anno sites need to be counted
    aac_sbs_key[cesa@mutations$amino_acid_change, pid := pid, on = 'aac_id']
    aac_sbs_key = aac_sbs_key[unique(multi_hits$pid), on = 'pid']
    
    if (cesa@maf[, .N] > 0) {
      maf_counts = cesa@maf[variant_type == 'sbs', .N, by = 'variant_id']
      aac_sbs_key[maf_counts, sbs_count := N, on = c(sbs_id = 'variant_id')]
    } else {
      aac_sbs_key[, sbs_count := 0]
    }
    aac_sbs_key[is.na(sbs_count), sbs_count := 0]
    maf_pid_counts = aac_sbs_key[, .(pid_freq = sum(sbs_count)), keyby = "pid"]
    multi_hits[maf_pid_counts, pid_freq := pid_freq, on = 'pid']
    multi_hits = merge.data.table(multi_hits, cesa@mutations$aac_sbs_key, by.x = 'variant_id', by.y = 'aac_id')
    multi_hits[, is_premature := aa_alt == "STOP" & aa_ref != "STOP"]
        
    # Any set of overlapping AACs has a single AAC chosen based on the following criteria:
    # MAF frequency (usually equal among all), essential splice status, premature stop codon, nonsilent status,
    # which protein has the most overall mutations in MAF data (will usually favor longer transcripts),
    # and finally just alphabetical on variant ID
    if(check_for_ref_data(cesa, 'transcripts')) {
        transcripts = get_ref_data(cesa, 'transcripts')
        multi_hits[transcripts, c('is_mane', 'is_mane_plus') := .(is_mane, is_mane_plus), on = c(pid = 'protein_id')]
        multi_hits = multi_hits[order(-essential_splice, -is_premature, aa_ref == aa_alt, -is_mane, -is_mane_plus, -maf_prevalence, -pid_freq, variant_id)]
        multi_hits[, c('is_mane', 'is_mane_plus') := NULL]
    } else {
        multi_hits = multi_hits[order(-essential_splice, -is_premature, aa_ref == aa_alt, -maf_prevalence, -pid_freq, variant_id)]
    }
    multi_hits[, is_premature := NULL]
    multi_hits[, sbs_id_dup := duplicated(sbs_id)]
    chosen_aac = multi_hits[, .(to_use = ! any(sbs_id_dup)), by = 'variant_id'][to_use == T, variant_id]
    
    # remove secondary (non-chosen) AACs, but save all sbs IDs and re-select those passing filters
    not_chosen_aac = setdiff(multi_hits$variant_id, chosen_aac)
    
    combined = combined[! not_chosen_aac, on = 'variant_id']
    
    sbs_in_chosen_aac = multi_hits[chosen_aac, sbs_id, on = 'variant_id'] # logically, no need for calling unique(sbs_id)
    sbs_in_not_chosen_aac = setdiff(multi_hits$sbs_id, sbs_in_chosen_aac)
    
    # We need to recover sbs that were in AACs and that passed user's filters,
    # but that are now no longer constituent sbs after de-overlapping.
    sbs_to_reselect = intersect(sbs_to_recover, sbs_in_not_chosen_aac)
    if (length(sbs_to_reselect) > 0) {
      reselected = select_variants(cesa, variant_ids = sbs_to_reselect)
      reselected = reselected[, .SD, .SDcols = names(combined)]
      combined = rbind(combined, reselected)
    }
  }
  # order output in chr/pos order
  setkey(combined, "chr")
  combined = setDT(combined[order(start)][BSgenome::seqnames(bsg), nomatch = NULL, on = "chr"])
  setcolorder(combined, c("variant_name", "variant_type", "chr", "start", "end", "variant_id", "ref", "alt", "gene", 
                          "strand", "aachange", "essential_splice", "intergenic", "trinuc_mut", "aa_ref", "aa_pos", "aa_alt", "coding_seq", 
                          "center_nt_pos", "pid", "multi_anno_site", 
                          "maf_prevalence", "samples_covering"))
  
  
  if(check_for_ref_data(cesa, "transcripts")) {
    transcripts = get_ref_data(cesa, "transcripts")
    combined[transcripts, let(is_MANE_transcript = is_mane,
                              transcript_tags = transcript_tags),
                              on = c(pid = 'protein_id')]
    combined[variant_type == 'sbs', is_MANE_transcript := NA]
  }

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






