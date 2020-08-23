#' Select and filter variants
#'
#' This is a lightweight function to help you search for and organization variant data
#' that is present in your CESAnalysis's MAF data and mutation annotation tables. You can
#' specify genes by name as well as CES-style variant IDs, and all variants matching
#' either will be returned. Or, leave both empty and all available variants will be
#' returned. The frequency of each variant in the MAF data will also be reported, and you
#' can filter output on this frequency.
#'
#' Note that while intergenic SNVs have their nearest genes annotated in the SNV tables,
#' these variants will not be captured by gene-based selection with this function, since
#' they're not actually in any gene.
#'
#' @param cesa CESAnalysis with MAF data loaded and annotated (e.g., with \code{load_maf()})
#' @param genes include variants (coding and noncoding) within these genes
#' @param variant_ids include variants
#' @param granges include all variants overlapping this GRanges object
#' @param min_freq filter out variants that occur in fewer than this many samples in the
#'   MAF data (Defaults to 1; set to 0 to include all passing variants that are in the
#'   annotation tables. Note that variants that are not in the annotation tables will
#'   never be returned. Use \code{generate_variants()} to annotate variants absent from the
#'   data.
#' @param include_subvariants Some mutations "contain" other mutations. For example, in
#'   cancereffectsizeR's, default hg19 reference data, KRAS_Q61H contains two constituent
#'   SNVs that both cause the same amino acid change: 12:25380275_T>G and 12:25380275_T>A.
#'   When include_subvariants = F (the default), and genes = "KRAS", output will be
#'   returned for KRAS_Q61H but not for the two SNVs (although their IDs will appear in
#'   the Q61H output). Set to true, and all three variants will be included in output, 
#'   assuming they don't get filtered out by min_freq or otherwise.
#' @param gr_padding add +/- this many bases to every range specified in "granges"
#'   (stopping at chromsome ends, naturally)
#' @export
select_variants = function(cesa, genes = NULL, variant_ids = NULL, granges = NULL, min_freq = 1, include_subvariants = F) {
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis object")
  }
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
  
  
  selected_snv_ids = character()
  selected_aac_ids = character()
  
  if (length(variant_ids) > 0) {
    variant_ids = unique(variant_ids)
    matching_snv_ids = cesa@mutations$snv[variant_ids, snv_id, nomatch = NULL]
    matching_aac_ids = cesa@mutations$amino_acid_change[variant_ids, aac_id, nomatch = NULL]
    missing_ids = setdiff(variant_ids, c(matching_snv_ids, matching_aac_ids))
    if (length(missing_ids) > 0) {
      num_missing = length(missing_ids)
      missing_ids = paste(missing_ids, collapse = ", ")
      stop(num_missing, " variants given in variant_ids couldn't be found. Either they're bad IDs or there are no annotations for them.\n",
           missing_ids)
    }
    selected_snv_ids = matching_snv_ids
    selected_aac_ids = matching_aac_ids
    
    # Add in SNV IDs that constitute the AACs
    if (include_subvariants) {
      selected_aac_ids = union(selected_aac_ids, cesa@mutations$amino_acid_change[matching_aac_ids, unlist(constituent_snvs)])
    }
  }
  
  
  if (length(genes) > 0) {
    genes = unique(genes)
    gene_names = get_genome_data(cesa, "gene_names")
    invalid_genes = genes[! genes %in% gene_names]
    num_invalid_genes = length(invalid_genes)
    if (num_invalid_genes > 0) {
      stop("Some of the selected genes do not appear in the CESAnalysis reference data:\n", paste(invalid_genes, collapse = ", "))
    }
    aac_in_genes = cesa@mutations$amino_acid_change[gene %in% genes, aac_id]
    
    # Note we're not returning intergenic SNVs that just have one of the chosen gene as their nearest gene
    genes_by_snv = cesa@mutations$snv[intergenic == FALSE, .(gene = unlist(genes)), by = "snv_id"]
    snv_in_genes = genes_by_snv[gene %in% genes, unique(snv_id)]
    
    selected_snv_ids = union(selected_snv_ids, snv_in_genes)
    selected_aac_ids = union(selected_aac_ids, aac_in_genes)
  }
  
  if (is.null(variant_ids) & is.null(genes)) {
    selected_aac_ids = cesa@mutations$amino_acid_change$aac_id
    selected_snv_ids = cesa@mutations$snv$snv_id
  }
  
  if (include_subvariants == FALSE) {
    subvariant_snv_ids = cesa@mutations$amino_acid_change[selected_aac_ids, unlist(constituent_snvs)]
    selected_snv_ids = setdiff(selected_snv_ids, subvariant_snv_ids)
  }
  
  selected_snv = cesa@mutations$snv[selected_snv_ids]
  selected_aac = cesa@mutations$amino_acid_change[selected_aac_ids]
  
  # Tabulate variants in MAF data and apply frequency filter
  snv_counts = cesa@maf[! is.na(snv_id), .N, by = "snv_id"][N >= min_freq]
  aac_counts = cesa@maf[! is.na(assoc_aac), .(aac_id = unlist(assoc_aac))][, .N, by = aac_id][N >= min_freq]
  selected_snv[, maf_frequency := 0]
  selected_snv = selected_snv[snv_counts, maf_frequency := N][maf_frequency >= min_freq]
  selected_snv[, assoc_aac2 := paste(unlist(assoc_aac), collapse = ","), by = "snv_id"]
  selected_snv[, assoc_aac := NULL][, assoc_aac := assoc_aac2][, assoc_aac2 := NULL]
  selected_aac[, maf_frequency := 0]
  selected_aac = selected_aac[aac_counts, maf_frequency := N][maf_frequency >= min_freq]
  
  
  # Create a combined, simplified table, with list columns collapsed
  selected_aac[, cs2 := paste(unlist(constituent_snvs), collapse = ","), by = "aac_id"]
  selected_aac[, constituent_snvs := NULL][, constituent_snvs := cs2][, cs2 := NULL]
  selected_aac[, variant_type := "AAC"]
  selected_aac[, intergenic := FALSE]
  selected_aac[, start := pmin(nt1_pos, nt3_pos)]
  selected_aac[, end := pmax(nt1_pos, nt3_pos)]
  selected_aac[, center_nt_pos := nt2_pos]
  selected_aac[, c("nt1_pos", "nt2_pos", "nt3_pos") := NULL]
  
  selected_snv[, variant_type := "SNV"]
  selected_snv[, constituent_snvs := NA_character_]
  selected_snv[, strand := NA_integer_] # because AAC table is +1/-1
  selected_snv[, c("start", "end") := .(pos, pos)]
  selected_snv[, pos := NULL]
  selected_aac[, assoc_aac := NA_character_]
  selected_aac[, all_genes := NA_character_]
  selected_aac[, multi_gene_hit := NA]
  
  selected_snv[, gene := sort(unlist(genes))[1], by = "snv_id"]
  
  selected_snv[, multi_gene_hit := F]
  is_multi_hit = sapply(selected_snv$genes, length) > 1
  selected_snv[, multi_gene_hit := is_multi_hit]
  notify_multi_genes = any(is_multi_hit)
  
  selected_snv[, all_genes := paste(unlist(genes), collapse = ","), by = "snv_id"]
  selected_snv[, genes := NULL]
  setnames(selected_snv, "snv_id", "variant_id")
  setnames(selected_aac, "aac_id", "variant_id")
  
  combined = rbindlist(list(selected_aac, selected_snv), use.names = T, fill = T)
  if(combined[, .N] == 0) {
    message("No variants passed selection criteria!")
    return(NULL)
  }
  
  cov_counts = cesa@samples[, .N, by = "covered_regions"]
  unique_combos = unique(combined$covered_in)
  setkey(cov_counts, "covered_regions")
  combo_counts = sapply(unique_combos, function(x) sum(cov_counts[x, N], na.rm = T))
  names(combo_counts) = sapply(unique_combos, function(x) paste0(x, collapse = "_"))
  combined[, maf_samples_covering := combo_counts[paste0(unlist(covered_in), collapse = "_")], by = "variant_id"]
  
  combined[, tmp := paste(unlist(covered_in), collapse = ","), by = "variant_id"]
  combined[, covered_in := NULL][, covered_in := tmp][, tmp := NULL]
  setkey(combined, "chr")
  
  # order output in chr/pos order
  bsg = BSgenome::getBSgenome(readRDS(paste0(cesa@genome_data_dir, '/genome_package_name.rds')))
  GenomeInfoDb::seqlevelsStyle(bsg) = "NCBI"
  
  combined = combined[BSgenome::seqnames(bsg), nomatch = NULL][, .SD[order(start)], by = "chr"]
  
  setcolorder(combined, c("variant_id", "variant_type", "chr", "start", "end", "ref", "alt", "gene", 
                          "strand", "aachange", "essential_splice", "intergenic", "assoc_aac", "trinuc_mut", "constituent_snvs",
                          "multi_gene_hit", "all_genes", "aa_ref", "aa_pos", "aa_alt", "coding_seq", "center_nt_pos", "pid", 
                          "covered_in", "maf_samples_covering", "maf_frequency"))
  
  
  if(notify_multi_genes) {
    message("Take note: Some of the returned SNVs have more than one gene/transcript annotation. There are comma-delimited\n",
            "in the \"all_genes\" column. The single gene in the \"gene\" field for these is just the first gene\n",
            "alphbetically (even if you selected variants by gene and didn't include this gene.) You can find these by\n",
            "filtering on the multi_gene_hit column.")
  }
  return(combined)
}

##' Add variant annotations
##' @param target_cesa CESAnalysis to receive variant annotations
##' @param snv_id character vector of CES-style snv_ids to validate, annotate (including
##'   finding and adding associated amino-acid-change mutations), and add to the target
##'   CESAnalysis
##' @param source_cesa source CESAnalysis whose snv_ids will be used to create new
##'   annotations in the target CESAnalysis
##' @export
# add_variants = function(target_cesa = NULL, snv_id = NULL, source_cesa = NULL) {
#   if(! is(target_cesa, "CESAnalysis")) {
#     stop("target_cesa should be a CESAnalysis", call. = F)
#   }
#   
#   if(! xor(is.null(snv_id, source_cesa))) {
#     stop("Add annotations via snv_id or source_cesa, but not both", call. = F)
#   }
#   
#   if (! is.null(source_cesa)) {
#     if(! is(source_cesa, "CESAnalysis")) {
#       stop("source_cesa should be a CESAnalysis", call. = F)
#     }
#     source_snv_table = source_cesa@mutations$snv
#     if (is.null(source_snv_table)) {
#       stop("source_cesa has no SNV annotations", call. = F)
#     }
#     snvs_to_annotate = source_cesa@mutations$snv$snv_id
#   }
#   
#   if(! is(snv_id, "character") | length(snv_id) == 0) {
#     stop("Expected snv_id to be character vector of snv_ids (e.g., 1:100_A>G", call. = F)
#   }
#   
#   # Load reference data if not already present
#   if (! target_cesa@ref_key %in% ls(.ces_ref_data)) {
#     preload_ref_data(target_cesa@ref_key)
#   }
#   
#   # will stop with errror if any IDs fail validation
#   validate_snv_ids(snv_ids, .ces_ref_data[[target_cesa@ref_key]]$genome)
#   
# }


#' validate_snv_ids
#' 
#' Ensures SNV IDs are valid for the given genome
#' 
#' @param snv_ids charcter vector of snv_ids
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

