#' Select and filter variants
#'
#' This function helps you find and view variant data from your CESAnalysis's MAF data and
#' mutation annotation tables. By default, almost all amino-acid-change mutations and
#' noncoding SNVs are returned. You can apply a series of filters to restrict output to
#' certain genes or genomic regions or require a minimum variant frequency in MAF data.
#' You can also specify some variants to include in output regardless of filters with
#' \code{variant_ids}. Special behavior: If \code{variant_ids} is used by
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
#'   \item intergenic: variant does not overlap any coding regions in the reference data
#'   \item trinuc_mut: for SNVs, the reference trinucleotide context, in deconstructSigs notation
#'   \item coding_seq: coding strand nucleotides in order of transcription
#'   \item center_nt_pos: regardless of strand, start/end give positions of two out of three AAC nucleotides; this
#'                       gives the position of the center nucleotide (maybe useful if the AAC spans a splice site)
#'   \item constituent_snvs: all SNVs that can produce a given variant 
#'   \item multi_anno_site: T/F whether variant has multiple gene/transcript/AAC annotations
#'   \item all_genes: all genes overlapping the variant in reference data
#'   \item covered_in: the names of all "covered regions" sets in the CESAnalysis that have coverage at the variant site
#'   \item maf_prevalence: number of occurrences of the variant in MAF data
#'   \item samples_covering: number of MAF samples with sequencing coverage at the variant site
#' }
#' 
#' @param cesa CESAnalysis with MAF data loaded and annotated (e.g., with \code{load_maf()})
#' @param genes Filter variants to specified genes.
#' @param min_freq Filter out variants with MAF frequency below threshold (default 0).
#'   Note that variants that are not in the annotation tables will never be returned. Use
#'   \code{add_variants()} to include variants absent from MAF data in your CESAnalysis.
#' @param variant_ids Vector of variant IDs to include in output regardless of
#'   filtering options. You can use CES-style AAC and SNV IDs or variant names like 
#'   "KRAS G12C". If this argument is used by itself (without any filtering arguments),
#'   then only these specified variants will be returned.
#' @param gr Filter out any variants not within input GRanges +/- \code{padding} bases.
#' @param variant_position_table Filter out any variants that don't intersect the
#'   positions given in chr/start/end of this table (1-based closed coordinates).
#'   Typically, the table comes from a previous \code{select_variants} call and can be expanded
#'   with \code{padding}. (Gritty detail: Amino acid change SNVs get special handling. Only the
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
#'   include_subvariants = T. Due to overlapping coding region definitions in reference
#'   data (e.g., genes with multiple transcripts), a site can have more than one
#'   amino-acid-change mutation annotation. To avoid returning the same genome-positional
#'   variants multiple times, the default is to return one AAC in these situations.
#'   Tiebreakers are MAF prevalence, essential splice site status, premature stop codon,
#'   non-silent status, gene/protein mutation count, alphabetical. If you set
#'   remove_secondary_aac to FALSE, you can't put the output variant table in selection
#'   calculation functions. An alternative is to set to FALSE, pick which
#'   (non-overlapping) variants you want, and then re-run select_variants() with those
#'   variants specified in \code{variant_ids}.
#' @return A data table with info on selected variants (see details), or a list of IDs.
#' @export
select_variants = function(cesa, genes = NULL, min_freq = 0, variant_ids = NULL, gr = NULL, variant_position_table = NULL, 
                            include_subvariants = F, padding = 0, collapse_lists = F, remove_secondary_aac = TRUE) {
  
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis object")
  }
  bsg = get_cesa_bsg(cesa)
  
  # Soon we'll check indels, too
  if(cesa@mutations$snv[, .N] == 0) {
    stop("There are no variants in the CESAnalysis.")
  }
  
  if(! is.numeric(min_freq) || length(min_freq) > 1 || min_freq - as.integer(min_freq) != 0) {
    stop("min_freq should be 1-length integer at least zero")
  }
  if(! is.null(genes) & ! is(genes, "character")) {
    stop("genes should be character vector of genes to include or NULL (for no gene filtering)")
  }
  
  if(! is.null(variant_ids) & ! is(variant_ids, "character")) {
    stop("variant_ids should be character vector of variant IDs to include, or left NULL.")
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

  # collect all variants, unless just variant_ids specified
  if (is.null(gr) && is.null(variant_position_table) && is.null(genes) && min_freq == 0 && length(variant_ids) > 0) {
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
    aac_passing_gr = cesa@mutations$aac_snv_key[captured$snv_id, unique(aac_id), on = 'snv_id', nomatch = NULL]
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
        msg = paste0("Shorthand amino-acid-change names (styled like \"KRAS_G12C\") were recognized and matched ",
                     "with cancereffectsizeR-style aac_ids. However, some of your variant names matched more than ",
                     "one aac_id (i.e., the same amino acid change is possible on multiple protein isoforms, and both are present ",
                     "in this analysis's variant annotations.) When the underlying SNVs are the same, only one AAC will be returned ",
                     "(unless you run with remove_secondary_aac = FALSE). Otherwise, all matching variants ",
                     "will be returned. (There may be additional matching variants with MAF frequency = 0 that are not ",
                     "annotated in this analysis; these will not be returned.)")
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
  
  
  # Rare situation: Constituent SNVs of AACs that get tossed in de-overlapping need to be
  # saved if they pass filters in their own right and don't overlap the AACs kept. Keep
  # track here of SNVS that need to appear in final output. Note if any of these don't
  # pass frequency filter, they shouldn't be saved, and they get removed below.
  snvs_to_recover = selected_snv_ids
  
  # Remove variants contained in other variants
  aac_snv_key = cesa@mutations$aac_snv_key[selected_aac_ids, on = 'aac_id']
  if (include_subvariants == FALSE) {
    selected_snv_ids = setdiff(selected_snv_ids, aac_snv_key$snv_id)
  }
  
  selected_snv = setDT(cesa@mutations$snv[selected_snv_ids])
  selected_aac = setDT(cesa@mutations$amino_acid_change[selected_aac_ids])
  
  # Get variant counts and coverage
  snv_from_aac = cesa@mutations$aac_snv_key[selected_aac$aac_id, .(aac_id, snv_id), on = 'aac_id']
  
  if(snv_from_aac[, .N] == 0 && length(selected_snv_ids) == 0) {
    message("No variants passed selection criteria!")
    return(NULL)
  }
  
  if (cesa@maf[, .N] > 0) {
    counts_and_cov = .variant_counts(cesa, samples = cesa@samples[, .(Unique_Patient_Identifier, covered_regions)],
                                     snv_from_aac = aac_snv_key[, .(aac_id, snv_id)],
                                     noncoding_snv_id = selected_snv_ids)
    setnames(counts_and_cov, c("total_prevalence", "total_covering"), c("maf_prevalence", "samples_covering"))
    
  } else {
    counts_and_cov = data.table(variant_id = c(selected_aac$aac_id, selected_snv_ids),
                                variant_type = c(rep('aac', selected_aac[, .N]), rep('snv', length(selected_snv_ids))),
                                maf_prevalence = 0, samples_covering = 0)
   
  }

  selected_snv[counts_and_cov, c("maf_prevalence", "samples_covering") := list(maf_prevalence, samples_covering), on = c(snv_id = 'variant_id')]
  selected_snv = selected_snv[maf_prevalence >= min_freq | snv_id %in% passlisted_ids]
  selected_aac[counts_and_cov, c("maf_prevalence", "samples_covering") := list(maf_prevalence, samples_covering), on = c(aac_id = 'variant_id')]
  selected_aac = selected_aac[maf_prevalence >= min_freq | aac_id %in% passlisted_ids]
  
  if (any(selected_snv$maf_prevalence < min_freq) || any(selected_aac$maf_prevalence < min_freq)) {
    pretty_message("Note: Some of your specifically-requested variants have MAF prevalence < min_freq. They will still appear in output.")
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
  selected_snv[cesa@mutations$aac_snv_key, multi_anno_site := multi_anno_site, on = 'snv_id']
  selected_snv[is.na(multi_anno_site), multi_anno_site := FALSE]
  setnames(selected_snv, c("genes", "snv_id"), c("all_genes", "variant_id"))
  
  # AACs get a short variant name that might not be uniquely identifying if a gene has more than one CDS
  selected_aac[, variant_name := paste(gene, aachange, sep = "_")]
  selected_aac[, variant_type := "aac"]
  selected_aac[, intergenic := FALSE]
  selected_aac[, start := pmin(nt1_pos, nt3_pos)]
  selected_aac[, end := pmax(nt1_pos, nt3_pos)]
  selected_aac[, center_nt_pos := nt2_pos]
  selected_aac[, c("nt1_pos", "nt2_pos", "nt3_pos") := NULL]
  
  if (selected_aac[, .N] > 0) {
    aac_to_snv = setDT(cesa@mutations$aac_snv_key[selected_aac$aac_id, on = 'aac_id'])
    # this form of assignment avoids data.table semantics and is faster on large list assignments
    aac_to_snv$genes = cesa@mutations$snv[aac_to_snv$snv_id, genes] 
    all_genes_by_aac_id = aac_to_snv[, .(genes = .(unique(unlist(genes)))), by = "aac_id"]
    selected_aac[all_genes_by_aac_id, all_genes := genes, on = "aac_id"]
    selected_aac[, nearest_pid := list(NA_character_)]
  } else {
    selected_aac[, all_genes := list(NA_character_)]
  }
  
  # We call an AAC a multi-anno site if any of its SNVs has multiple annotations.
  selected_aac[cesa@mutations$aac_snv_key, multi_anno_site := any(multi_anno_site), on = 'aac_id']
  setnames(selected_aac, "aac_id", "variant_id")

  # Combine SNV and AAC tables
  combined = rbindlist(list(selected_aac, selected_snv), use.names = T, fill = T)
  if(combined[, .N] == 0) {
    message("No variants passed selection criteria!")
    return(NULL)
  }
  # convert 0-length covered_in to NA (for sites just covered in whole-genome)
  combined[which(sapply(covered_in, length) == 0), covered_in := list(NA_character_)]
  
  # collapse list columns, if specified
  if (collapse_lists) {
    # Problem: unstrsplit converts NA to "NA"
    # "NA" is not a valid value for any of these except all_genes, and going to assume there is not a gene called "NA"
    list_cols = c("constituent_snvs", "covered_in", "all_genes")
    combined[, (list_cols) := lapply(.SD, function(x) S4Vectors::unstrsplit(x, sep = ",")), .SDcols = list_cols]
    combined[, (list_cols) := lapply(.SD, function(x) gsub('NA', NA_character_, x)), .SDcols = list_cols]

  }
  
  # handle overlapping  mutations using tiebreakers explained below
  if (remove_secondary_aac) {
    multi_hits = combined[variant_type == "aac" & multi_anno_site == TRUE]
    num_to_check = multi_hits[, .N]
    if (num_to_check > 0) {
      # for tie-breaking, count how many mutations are in each gene found in these multi_hit records
      # will need to produce these counts from scratch since some of the variants may not be in this select_variants() run
      aac_snv_key = cesa@mutations$aac_snv_key[multi_anno_site == TRUE] # only multi-anno sites need to be counted
      aac_snv_key[cesa@mutations$amino_acid_change, pid := pid, on = 'aac_id']
      aac_snv_key = aac_snv_key[unique(multi_hits$pid), on = 'pid']
      
      if (cesa@maf[, .N] > 0) {
        maf_counts = cesa@maf[variant_type == 'snv', .N, by = 'variant_id']
        aac_snv_key[maf_counts, snv_count := N, on = c(snv_id = 'variant_id')]
      } else {
        aac_snv_key[, snv_count := 0]
      }
      aac_snv_key[is.na(snv_count), snv_count := 0]
      maf_pid_counts = aac_snv_key[, .(pid_freq = sum(snv_count)), keyby = "pid"]
      multi_hits[maf_pid_counts, pid_freq := pid_freq, on = 'pid']
      
      # Any set of overlapping AACs has a single AAC chosen based on the following criteria:
      # MAF frequency (usually equal among all), essential splice status, premature stop codon, nonsilent status,
      # which protein has the most overall mutations in MAF data (will usually favor longer transcripts),
      # and finally just alphabetical on variant ID
      multi_hits[, is_premature := aa_alt == "STOP" & aa_ref != "STOP"]
      multi_hits = multi_hits[order(-maf_prevalence, -essential_splice, -is_premature, aa_ref == aa_alt, -pid_freq, variant_id)]
      multi_hits[, is_premature := NULL]
      
      # When chr/nt/aachange all match, the higher-up entry in the table will always be
      # chosen to the exclusion of other matches
      original_multi_hit_ids = multi_hits$variant_id
      original_const_snv = multi_hits$constituent_snvs
      multi_hits = unique(multi_hits, by = c("chr", "start", "end", "center_nt_pos", "aa_alt"))
      setkey(multi_hits, 'variant_id', physical = F) # important not to re-sort since we just sorted
      chosen_aac = new.env(parent = emptyenv(), size = 30 * multi_hits[, .N])
      covered_snv = new.env(parent = emptyenv(), size = 60 * multi_hits[, .N])
      mapply(
        function(curr_candidate, curr_covered_snv) {
          # We will use the current AAC only if none of the constituent SNVs have been used yet
          for (i in curr_covered_snv) {
            if(exists(i, covered_snv)) {
              return()
            }
          }
          for (i in curr_covered_snv) {
            covered_snv[[i]] = TRUE
          }
          chosen_aac[[curr_candidate]] = TRUE
        }, multi_hits$variant_id, multi_hits$constituent_snvs)
      
      # remove secondary (non-chosen) AACs, but save all SNV IDs and re-select those passing filters
      chosen_aac = ls(chosen_aac, sorted = FALSE)
      not_chosen_aac = setdiff(original_multi_hit_ids, chosen_aac)
      remaining_const_snv = multi_hits[chosen_aac, unlist(constituent_snvs), on = 'variant_id']
      combined = setDT(combined[! not_chosen_aac, on = 'variant_id'])
      
      # Edge case: Will need to get annotations for SNVs that are in AAC, that passed user's filters,
      # but that are now no longer constituent SNVs after de-overlapping
      snv_to_reselect = intersect(snvs_to_recover, setdiff(original_const_snv, remaining_const_snv))
      if (length(snv_to_reselect) > 0) {
        reselected = select_variants(cesa, variant_ids = snv_to_reselect)[maf_prevalence >= min_freq]
        combined = rbind(combined, reselected)
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
  combined = setDT(combined[order(start)][BSgenome::seqnames(bsg), nomatch = NULL, on = "chr"])
  setcolorder(combined, c("variant_name", "variant_type", "chr", "start", "end", "variant_id", "ref", "alt", "gene", 
                          "strand", "aachange", "essential_splice", "intergenic", "nearest_pid", "trinuc_mut", "aa_ref", "aa_pos", "aa_alt", "coding_seq", 
                          "center_nt_pos", "pid", "constituent_snvs", "multi_anno_site", "all_genes",
                          "covered_in", "maf_prevalence", "samples_covering"))
  
  setattr(combined, "cesa_id", cesa@advanced$uid)
  setattr(combined, "nonoverlapping", nonoverlapping)
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






