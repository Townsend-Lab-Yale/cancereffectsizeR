#' Add variant annotations
#' 
#' Use this function to add variant annotations to your CESAnalysis by specifying variants
#' to add in one of five ways: a data.table containing genomic coordinates (output from
#' select_variants(), typically), a GRanges object, a BED file, another CESAnalysis, or
#' SNV/AAC IDs.
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
#' @param aac_id Character vector of AAC IDs (or short names, like "KRAS_G12C")
#' @param source_cesa Another CESAnalysis from which to copy snv_ids. SNVs will be
#'   re-annotated using the target_cesa's associated reference data.
#' @param padding How many bases (default 0) to expand start and end of each gr range
#' @export
add_variants = function(target_cesa = NULL, variant_table = NULL, snv_id = NULL, aac_id = NULL, bed = NULL, 
                        gr = NULL, source_cesa = NULL, padding = 0) {
  if(! is(target_cesa, "CESAnalysis")) {
    stop("target_cesa should be a CESAnalysis", call. = F)
  }
  target_cesa = copy_cesa(target_cesa)
  target_cesa = update_cesa_history(target_cesa, match.call())
  
  if (! is(padding, "numeric") || length(padding) != 1 || trunc(padding) != padding || padding < 0) {
    stop("padding should be 1-length non-negative integer")
  }
  
  # just one possible method should be chosen
  if (sum(sapply(list(variant_table, gr, bed, snv_id, aac_id, source_cesa), is.null)) != 5){
    stop("Use exactly one option out of variant_table, gr, bed, snv_id, aac_id, source_cesa.")
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
      if(sum(sapply(target_cesa@mutations, nrow)) == 0) {
        target_cesa@mutations = copy(source_cesa@mutations)
      } else {
        new_snvs = source_cesa@mutations$snv[! target_cesa@mutations$snv$snv_id, on = "snv_id"]
        new_dbs = source_cesa@mutations$dbs[! target_cesa@mutations$dbs$dbs_id, on = "dbs_id"]
        if (new_snvs[, .N] + new_dbs[, .N] == 0) {
          warning("There were no new variants to copy over.")
        }
        # covered_in may vary, but doesn't matter because it will be regenerated from scratch
        target_cesa@mutations$snv = rbind(target_cesa@mutations$snv, new_snvs)
        target_cesa@mutations$dbs = rbind(target_cesa@mutations$dbs, new_dbs)
        
        if(source_cesa@mutations$amino_acid_change[, .N] > 0) {
          new_aacs = source_cesa@mutations$amino_acid_change[! target_cesa@mutations$amino_acid_change$aac_id, on = "aac_id"]
          target_cesa@mutations$amino_acid_change = rbind(target_cesa@mutations$amino_acid_change, new_aacs)
          target_cesa@mutations$aac_snv_key = unique(rbind(target_cesa@mutations$aac_snv_key,
                                                           source_cesa@mutations$aac_snv_key))
        }
        if (new_dbs[, .N] > 0) {
          target_cesa@mutations$dbs = unique(rbind(target_cesa@mutations$dbs,
                                                           source_cesa@mutations$dbs), by = 'dbs_id')
        }
        if(source_cesa@mutations$dbs_codon_change[, .N] > 0) {
          target_cesa@mutations$dbs_codon_change = unique(rbind(target_cesa@mutations$dbs_codon_change,
                                                           source_cesa@mutations$dbs_codon_change), by = 'dbs_aac_id')
        }
        target_cesa@mutations$aac_dbs_key = unique(rbind(target_cesa@mutations$aac_dbs_key,
                                                         source_cesa@mutations$aac_dbs_key))
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
      if(! any(ref %in% nt)) {
        pretty_message("All new variants are at sites with ambiguous reference sequence, so no annotations can be added.")
        return(target_cesa)
      }
      snv_table = snv_table[ref %in% nt]
      message("Note: Some variants in input were dropped because of ambiguous reference sequence (N's)")
    }
    snv_table = snv_table[rep(1:.N, each = 4)]
    snv_table[, alt := rep.int(c("A", "C", "G", "T"), .N/4)]
    snv_table = snv_table[ref != alt]
    snvs_to_annotate = snv_table[, paste0(chr, ':', pos, '_', ref, '>', alt)]
  }
  
  # If supplied SNV IDs (rather than source_cesa, gr, bed, variant_table), validate them
  if(! is.null(snv_id)) {
    if(! is(snv_id, "character") || length(snv_id) == 0) {
      stop("Expected snv_id to be character vector of snv_ids (e.g., 1:100_A>G", call. = F)
    }
    # will stop with errror if any IDs fail validation
    validate_snv_ids(snv_id, get_cesa_bsg(target_cesa))
    snvs_to_annotate = snv_id
  }
  
  # If supplied aac_id, get constituent SNVs
  refset = .ces_ref_data[[target_cesa@ref_key]]
  if (! is.null(aac_id)) {
    if(! is.character(aac_id) || length(aac_id) == 0) {
      stop("Expected aac_id to be a character vector of amion-acid-change IDs (or variant names like KRAS_G12C).")
    }
    aac_id = complete_aac_ids(partial_ids = aac_id, refset = refset)
    
    if(length(aac_id) == 0) {
      stop('No valid aac_id to add.')
    }
    aac_dt = setDT(tstrsplit(aac_id, split = "_"))
    setnames(aac_dt, c('gene', 'aachange', 'pid'))
    if ('gene' %in% names(GenomicRanges::mcols(refset$gr_genes))) {
      aac_dt[, entry_name := pid]
    } else {
      aac_dt[, entry_name := gene]
    }
    
    aac_dt[, char_length := nchar(aachange)]
    aac_dt[, aa_alt := substr(aachange, char_length, char_length)]
    aac_dt[, aa_pos := as.numeric(substr(aachange, 2, char_length - 1))]
    
    # for compatibility with internal codon table, convert to three-letter codes
    aac_dt[, aa_alt := seqinr::aaa(aa_alt)]
    aac_dt[aa_alt == 'Stp', aa_alt := 'STOP'] # also for compatibility
    
    fast_refcds = list2env(refset$RefCDS[unique(aac_dt$entry_name)])
    
    
    snvs_to_annotate = unique(unlist(mapply(aac_to_snv_ids, 
                                     aa_pos = aac_dt$aa_pos,
                                     aa_alt = aac_dt$aa_alt,
                                     refcds_entry_name = aac_dt$entry_name,
                                     MoreArgs = list(bsg = refset$genome, refcds = fast_refcds), SIMPLIFY = F)))
  }
  
  num_variants = length(snvs_to_annotate)
  if(num_variants == 0) {
    stop("No variants to add (check your input).")
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
  
  # split snv_ids into MAF-like table for annotation
  cesa = target_cesa
  snv_id = snvs_to_annotate
  maf = as.data.table(tstrsplit(snv_id, split = '[:_>]'))
  colnames(maf) = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")
  maf[, Start_Position := as.numeric(Start_Position)]
  
  annotations = annotate_variants(refset = .ces_ref_data[[target_cesa@ref_key]], variants = maf)
  aac_table = annotations$amino_acid_change
  snv_table = annotations$snv
  if (aac_table[, .N] > 0) {
    cesa@mutations[["amino_acid_change"]] = unique(rbind(cesa@mutations$amino_acid_change, aac_table, fill = T), by = "aac_id")
    setkey(cesa@mutations$amino_acid_change, "aac_id")
  }
  cesa@mutations[["snv"]] = unique(rbind(cesa@mutations$snv, snv_table, fill = T), by = "snv_id")
  cesa@mutations[["aac_snv_key"]] = unique(rbind(cesa@mutations$aac_snv_key, annotations$aac_snv_key))
  cesa@mutations[["dbs"]] = unique(rbind(cesa@mutations$dbs, annotations$dbs))
  cesa@mutations[["dbs_codon_change"]] = unique(rbind(cesa@mutations$dbs_codon_change, annotations$dbs_codon_change))
  setkey(cesa@mutations$aac_snv_key, 'aac_id')
  setkey(cesa@mutations$snv, "snv_id")
  
  # Record which coverage ranges each new variant is covered in
  cesa@advanced$add_variants_used = TRUE # if add_variants has run, can't take shortcuts in update_covered_in 
  cesa = update_covered_in(cesa)
  return(cesa)
}

#' Get SNVs that cause an amino acid change
#' 
#' An internal function to figure out the SNVs that can cause a given amino acid substitution in a transcript
#' 
#' @param refcds_entry A RefCDS entry for the relevant transcript
#' @param aa_pos Integer position of substitution on the transcript.
#' @param aa_alt Identity of substitution, either a three-letter code ("Lys") or "STOP"
#' @param bsg A BSgenome object for the genome build associated with the RefCDS entry
#' @keywords internal
aac_to_snv_ids = function(refcds_entry_name, aa_pos, aa_alt, bsg, refcds) {
  entry = refcds[[refcds_entry_name]]
  chr = entry$chr
  neg_strand = entry$strand == -1
  
  if (neg_strand) {
    strand = '-'
    cds_int = as.data.table(entry$intervals_cds)[.N:1][, .(start = V2, end = V1)]
  } else {
    strand = '+'
    cds_int = as.data.table(entry$intervals_cds)[, .(start = V1, end = V2)]
  }
  
  # From the amino acid position, get the genomic positions of the three nucleotides.
  transcript_pos = 3*(aa_pos - 1) + c(1, 2, 3)
  genome_pos = unlist(mapply(seq, cds_int$start, cds_int$end))[transcript_pos]
  
  # Get the codon and positive-strand reference sequence at the positions 
  codon = unlist(getSeq(bsg, GPos(seqnames = chr, pos = genome_pos, strand = strand)))
  ref_seq = codon
  if (neg_strand) {
    ref_seq = Biostrings::complement(codon)
  }
  
  to_annotate = codon_snvs_to_aa[[as.character(codon)]][[aa_alt]]
  snv_id = character()
  for (i in 1:3) {
    if(length(to_annotate[[i]]) > 0) {
      dna_alt =to_annotate[[i]]
      if(neg_strand) {
        dna_alt = seqinr::comp(dna_alt, forceToLower = F)
      }
      snv_id = c(snv_id, paste0(chr, ':', genome_pos[i], '_', as.character(ref_seq[i]), '>', dna_alt))
    }
  }
  return(snv_id)
}
