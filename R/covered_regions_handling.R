
#' add_covered_regions
#' @param target_cesa CESAnalysis with annotated variants that the covered regions will be added to
#' @param source_cesa Another CESAnalysis to copy all covered regions from
#' @param covered_regions A GRanges object or BED file path with genome build matching the target_cesa,
#'                      if not using source_cesa
#' @param covered_regions_name A name to identify the covered regions, if not using source_cesa
#' @param coverage_type exome, genome, or targeted (if not using source_cesa)
#' @param covered_regions_padding optionally, add +/- this many bp to each interval in covered_regions
#' @return CESAnalysis given in target_cesa, with the new covered regions added
#' @export
add_covered_regions = function(target_cesa = NULL, source_cesa = NULL, covered_regions = NULL, 
                               covered_regions_name = NULL, coverage_type = NULL, covered_regions_padding = 0) {

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
  
  target_cesa = copy_cesa(target_cesa)
  prev_recording_status = target_cesa@advanced$recording
  target_cesa = update_cesa_history(target_cesa, match.call())
  target_cesa@advanced$recording = F
  
  # Validate both possible sets of arguments
  if (! is.null(source_cesa)) {
    if (! is(source_cesa, "CESAnalysis")) {
      stop("source_cesa should be a CESAnalysis")
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
                                          covered_regions_name = tgs_set, coverage_type = "targeted")
    }
    
    wgs_sets = names(source_cesa@coverage$genome)
    for (wgs_set in wgs_sets) {
      target_cesa = assign_gr_to_coverage(target_cesa, gr = source_cesa@coverage$genome[[tgs_set]],
                                          covered_regions_name = wgs_set, coverage_type = "genome")
    }
    
    
    target_cesa@advanced$recording = prev_recording_status
    return(update_covered_in(target_cesa))
  } else {
    target_cesa@advanced$recording = prev_recording_status
    return(.add_covered_regions(cesa = target_cesa, coverage_type = coverage_type, covered_regions = covered_regions,
                                covered_regions_name = covered_regions_name,covered_regions_padding = covered_regions_padding))
  }
}

#' .add_covered_regions
#' 
#' @param covered_regions A GRanges object or BED file path with genome build matching the target_cesa,
#'                      if not using source_cesa
#' @param covered_regions_name A name to identify the covered regions, if not using source_cesa
#' @param coverage_type exome or targeted, if not using source_cesa
#' @param covered_regions_padding optionally, add +/- this many bp to each interval in covered_regions
#' @keywords internal
#' @return CESAnalysis given in target_cesa, with the new covered regions added
.add_covered_regions = function(cesa, coverage_type, covered_regions, covered_regions_name, covered_regions_padding) {
  if (! is.character(coverage_type) | length(coverage_type) != 1 | ! coverage_type %in% c("exome", "genome", "targeted")) {
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
  legal_name = '^[a-z][-0-9a-z\\_\\.]*$'
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
  gr = clean_granges_for_cesa(cesa = cesa, gr = gr, padding = covered_regions_padding)
  cesa = assign_gr_to_coverage(cesa, gr = gr, covered_regions_name = covered_regions_name, coverage_type = coverage_type)
  cesa = update_covered_in(cesa)
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
  if (coverage_type %in% c("exome", "genome") && check_for_ref_data(cesa, "default_exome_gr")) {
    covered_regions_bases_covered = sum(IRanges::width(IRanges::ranges(gr)))
    generic_bases_covered = sum(IRanges::width(IRanges::ranges(get_ref_data(cesa, "default_exome_gr"))))
    if (coverage_type == "exome") {
      if (covered_regions_bases_covered / generic_bases_covered < .4) {
        warning(paste0("Input coverage intervals are described as whole-exome but are less than 40% of the size of this genome's default exome intervals.\n",
                       "This might make sense if your exome capture array is  lean, but if this is actually targeted sequencing data,\n",
                       "start over with the coverage=\"targeted\"."))
      }
    } else {
      # taking genome size to be the sum of the lengths of all CES-supported contigs
      genome_size = sum(GenomeInfoDb::seqlengths(get_cesa_bsg(cesa))[cesa@advanced$genome_info$supported_chr])
      if (covered_regions_bases_covered / genome_size < .4) {
        msg = paste0("Input coverage intervals are described whole-genome but cover less than 40% of the genome. This might make ",
               "sense if you are excluding genomic regions due to repetitiveness or mappability, although excluding ",
               "consensus repetitive regions on the human genome would not be this restrictive. Consider whether the coverage regions ",
               "are better described as exome or targeted.")
        msg = pretty_message(msg, emit = F)
        warning(msg)
      }
    }
  }
  
  if (! coverage_type %in% names(cesa@coverage)) {
    cesa@coverage[[coverage_type]] = list()
  }
  cesa@coverage[[coverage_type]][[covered_regions_name]] = gr
  return(cesa)
}


#' update_covered_in
#' 
#' Updates the covered_in annotation for all variants
#' to include all covered regions in the CESAnalysis
#' 
#' Also updates internal cached output of select_variants().
#' 
#' @param cesa CESAnalysis
#' @return CESAnalysis with regenerated covered-in annotations
#' @keywords internal
update_covered_in = function(cesa) {
  # if no SNVs, nothing to update (may change once indels are supported)
  if(cesa@mutations$snv[, .N] == 0) {
    return(cesa)
  }
  
  if(length(cesa@coverage) == 0) {
    return(cesa)
  }
  
  snv_table = cesa@mutations$snv
  aac_table = cesa@mutations$amino_acid_change
  

  
  # simpler to delete old annotations and re-do rather than updating
  snv_table[, covered_in := NULL]

  
  snv_gr = GenomicRanges::makeGRangesFromDataFrame(snv_table, seqnames.field = "chr", start.field = "pos", end.field = "pos")
  
  # test each MAF locus against all coverage grs
  # this returns a data frame where rows match MAF rows, columns are T/F for each coverage gr
  all_coverage = c(cesa@coverage[["exome"]], cesa@coverage[["targeted"]], cesa@coverage[["genome"]]) # names shouldn't overlap due to add_covered_regions validation
  
  # When exome+ is being used, not much point of including the original "exome" ranges
  if ('exome+' %in% names(all_coverage)) {
    all_coverage = all_coverage[names(all_coverage) != 'exome']
  }
  
  # all_coverage will be NULL if there's only full-coverage WGS data
  just_exome_plus = FALSE
  if(is.null(all_coverage)) {
    snv_table[, covered_in := list()]
  } else if (length(all_coverage) == 1 && names(all_coverage) == 'exome+' && 
             sum(cesa@samples$covered_regions == 'genome') == 0 && is.null(cesa@advanced$add_variants_used)) {
    # If only exome+ data and add_variants() hasn't been used to import (potentially uncovered) variants,
    # no need to do any calculations.
    #
    # Oddity: Consider out-of-exome AACs with multiple constituent SNVs but not all
    # present in MAF data. Those SNVs get discovered in annotate_variants(), after exome+
    # coverage has been defined. However, this shortcut will update them to be covered.
    # But if a user has used add_variants() previously, those SNVs will stay uncovered.
    # However, this isn't much of an issue because the AACs are still considered covered,
    # and the issue only occurs on fringe AACs that aren't part of default exome coverage.
    snv_table[, covered_in := list(list('exome+'))]
    just_exome_plus = TRUE
  } else {
    all_coverage = all_coverage[sort(names(all_coverage))] # sort by covered regions name
    is_covered = setDT(lapply(all_coverage, function(x) IRanges::overlapsAny(snv_gr, x, type = "within")))
    all_names = names(is_covered)
    # when we support just R > 4.1, can just use simplify = F on apply to avoid playing these games
    covered_in = apply(is_covered, 1, function(x) all_names[x])
    if(is(covered_in, "character")) {
      covered_in = as.list(covered_in)
    } else if(is(covered_in, "matrix")) {
      covered_in = lapply(1:ncol(covered_in), function(x) covered_in[, x])
    }
    snv_table[, covered_in := list(covered_in)]
  }
  cesa@mutations$snv = snv_table
  setkey(cesa@mutations$snv, 'snv_id')
  
  # We're going to cheat and say samples have coverage at aac sites if they have coverage on any of the three codon positions
  # in practice, there probably is coverage even if the coverage intervals disagree
  if (! is.null(aac_table)) {
    if (just_exome_plus) {
      aac_table[, covered_in := list(list('exome+'))]
    } else {
      if ("covered_in" %in% colnames(aac_table)) {
        aac_table[, covered_in := NULL]
      }
      aac_coverage = cesa@mutations$aac_snv_key[cesa@mutations$amino_acid_change$aac_id, on = 'aac_id']
      aac_coverage[snv_table, covered_in := covered_in, on = 'snv_id']
      aac_coverage = aac_coverage[, .(covered_in = list(unique(unlist(covered_in)))), by = 'aac_id']
      
      # edge case of adding single variant to empty CESAnalysis
      if(aac_coverage[, .N] == 1 && is.null(aac_coverage$covered_in[[1]])) {
        aac_table[, covered_in := list(NA_character_)]
      } else if(aac_coverage[, .N] > 0) {
        aac_table[aac_coverage, covered_in := covered_in, on = "aac_id"]
      }
    }
    cesa@mutations$amino_acid_change = aac_table
    setkey(cesa@mutations$amino_acid_change, 'aac_id')
  }
  
  # Finally, update internal cached variants
  cesa@advanced$cached_variants = suppressMessages(select_variants(cesa))
  return(cesa)
}
