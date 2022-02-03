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

#' Create full AAC ID
#' 
#' For example, KRAS_G12C -> KRAS_G12C_ENSP00000256078 (ces.refset.hg19). In cases of
#' multiple protein IDs per gene, will return more IDs than input. Otherwise,
#' input/output will maintain order.
#' 
#' Output IDs are not fully validated. Use validate_aac_ids().
#' 
#' @param partial_ids AAC variant names, such as "KRAS_G12C" or "MIB2 G395C"
#' @param refset reference data set (environment object)
#' @keywords internal
complete_aac_ids = function(partial_ids, refset) {
  if(! is.environment(refset) || is.null(refset[["gr_genes"]])) {
    stop("refset doesn't appear to be a CES reference data set")
  }
  
  if (! is.character(partial_ids)) {
    stop("partial_ids expected to be character")
  }
  
  if(anyNA(partial_ids)) {
    stop("There are NA entries in partial_ids")
  }
  
  # Not carefully validating, but if the IDs all look complete already, return them
  looks_complete = grepl('.+_.+_', partial_ids)
  if (all(looks_complete)) {
    problems = validate_aac_ids(partial_ids, refset)
    if(length(problems) > 0) {
      print(problems)
      stop('Input IDs look like full AAC IDs (gene_aachange_pid), but some are invalid (see above).')
    }
    return(partial_ids)
  }
  if (any(looks_complete)) {
    stop('Encountered an apparent mix of partial and full AAC IDs.')
  }
  
  # replace a space to allow things like "BRAF V600E"
  gene_and_aachange = gsub(" ", "_", partial_ids)
  dt = data.table(gene_and_aachange = gene_and_aachange)
  dt[, gene := sub('_.*', '', gene_and_aachange)]
  
  invalid_genes = setdiff(dt$gene, refset$gene_names)
  if(length(invalid_genes) > 0) {
    stop("Invalid genes (expected gene_aachange):\n", paste(invalid_genes, collapse = ", "))
  }
  
  refcds = refset$RefCDS
  gr_cds = refset$gr_genes
  if ("gene" %in% names(GenomicRanges::mcols(gr_cds))) {
    gene_to_pid = unique(as.data.table(GenomicRanges::mcols(gr_cds)))
    setnames(gene_to_pid, 'names', 'protein_id')
  } else {
    gene_to_pid = rbindlist(lapply(refcds[unique(dt$gene)], '[', 'protein_id'), idcol = 'gene')
  }
  dt = gene_to_pid[dt, on = 'gene'] # take all combinations of pid and gene_aachange (even though some will be invalid)
  tentative_ids = dt[, paste(gene_and_aachange, protein_id, sep = '_')]
  # determine which AAC are valid (verifies ref/alt amino acids are correct and possible at given position)
  bad_ids = unlist(validate_aac_ids(aac_ids = tentative_ids, refset = refset))
  good_ids = setdiff(tentative_ids, bad_ids)
  
  final_gene_and_aachange = unique(gsub('(.*)_[^_]+$', '\\1', good_ids))
  not_completed = setdiff(gene_and_aachange, final_gene_and_aachange)
  if(length(not_completed) > 0) {
    warning(paste0('Some input AAC variant names could not be expanded to valid AAC IDs, so no returned IDs correspond to them:\n',
                          paste(not_completed, sep = ', '), '.'))
  }
  return(good_ids)
}


#' Ensure AAC IDs are valid for a given reference data set
#' 
#' Given a vector of AAC IDs, determines whether each is valid/possible, returning
#' NULL if all are valid, or a list of problems, or an error if input isn't parseable.
#' 
#' An ID is invalid if the gene and/or protein_id are not in the reference data, or if the
#' reference amino acid ("G" in KRAS_G12C) is incorrect, or if there is no possible SNV
#' that can create the proposed change. For example, KRAS G12C is possible when the codon
#' (GGT) acquires a G>T substitution in the first position, while G12K is not possible
#' because no single substitution can transform GGT into a lysine codon.
#' 
#' @param aac_ids AAC variant IDs
#' @param refset reference data set (environment object)
#' @keywords internal
validate_aac_ids = function(aac_ids, refset) {
  if(! is.environment(refset) || is.null(refset[["gr_genes"]])) {
    stop("refset doesn't appear to be a CES reference data set")
  }
  
  if (! is.character(aac_ids)) {
    stop("aac_ids expected to be character")
  }
  if(anyNA(aac_ids)) {
    stop("aac_ids contains NA values.")
  }
  
  if (length(aac_ids) == 0) {
    return(NULL) # No IDs = no problem
  }
  dt = setDT(tstrsplit(aac_ids, '_'))
  
  
  if(length(dt) != 3 || anyNA(dt)) {
    stop('One or more misformatted AAC IDs: expected gene_aachange_pid.')
  }
  
  setnames(dt, c("gene", "aachange", "pid"))
  dt[, input_id := aac_ids]
  
  
  problems = list()
  bad_gene = dt[! gene %in% refset[["gene_names"]], which = T]
  
  if(length(bad_gene) > 0) {
    problems[["no_such_gene"]] = dt[bad_gene, input_id]
    dt = dt[! bad_gene]
  }
  
  
  gr_cds = refset[["gr_genes"]]
  refcds = refset[["RefCDS"]]
  if ("gene" %in% names(GenomicRanges::mcols(gr_cds))) {
    dt[, entry_name := pid]
    refcds = refcds[unique(dt$entry_name)]
    cds_info = rbindlist(lapply(refcds, '[', 
                                c("real_gene_name", "CDS_length")), idcol = 'entry_name')
    dt[cds_info, expected_gene := real_gene_name, on = 'entry_name']
    pid_gene_mismatch = dt[expected_gene != gene, which = T]
  } else {
    dt[, entry_name := gene]
    refcds = refcds[unique(dt$entry_name)]
    cds_info = rbindlist(lapply(refcds[unique(dt$entry_name)], '[', 
                                c("protein_id", "CDS_length")), idcol = 'entry_name')
    dt[cds_info, expected_pid := protein_id, on = 'entry_name']
    pid_gene_mismatch = dt[expected_pid != pid, which = T]
  }
  
  if(length(pid_gene_mismatch) > 0) {
    problems[["gene_pid_mismatch"]] = dt[pid_gene_mismatch, input_id]
    dt = dt[! pid_gene_mismatch]
  }
  
  cds_info[, seq_cds := list(sapply(refcds, '[[', 'seq_cds'))]
  
  dt[, tmp_length := nchar(aachange)]
  dt[, aa_ref := substr(aachange, 1, 1)]
  dt[, aa_alt := substr(aachange, tmp_length, tmp_length)]
  dt[, aa_pos := substr(aachange, 2, tmp_length - 1)]
  
  dt[, aa_pos := suppressWarnings(as.integer(aa_pos))]
  noninteger_aa_pos = dt[is.na(aa_pos), which = T]
  if(length(noninteger_aa_pos) > 0) {
    problems[["non-integer_aa_pos"]] = dt[noninteger_aa_pos, input_id]
    dt = dt[! noninteger_aa_pos]
  }
  
  dt[, seq_start := 3 * aa_pos - 2]
  dt[cds_info, c("cds_length", "seq_cds") := list(CDS_length, seq_cds), on = "entry_name"]
  pos_oob = dt[cds_length < seq_start, which = T]
  if(length(pos_oob) > 0) {
    problems[["aa_pos_out-of-bounds"]] = dt[pos_oob, input_id]
    dt = dt[! pos_oob]
  }
  codons = Biostrings::subseq(DNAStringSet(dt$seq_cds), start = dt$seq_start, width = 3)
  dt[, actual_ref := as.character(Biostrings::translate(codons, no.init.codon = T))]
  incorrect_ref = dt[aa_ref != actual_ref, which = T]
  if(length(incorrect_ref) > 0) {
    problems[['incorrect_aa_ref']] = dt[incorrect_ref, input_id]
    dt = dt[! incorrect_ref]
    codons = codons[-incorrect_ref] # ! doesn't work
  }
  
  # drop synonymous records (allowing those for now)
  is_synonymous = dt[aa_ref == aa_alt, which = T]
  if(length(is_synonymous) > 0) {
    dt = dt[! is_synonymous]
    codons = codons[-is_synonymous]
  }
  
  dt[aa_alt %in% c('*', 'STP'), aa_alt := 'STOP']
  dt[aa_alt != 'STOP', aa_alt := suppressWarnings(seqinr::aaa(aa_alt))] # for compatibility with internal table
  
  # handle aa_alt that are not valid amino acid abbreviations
  invalid_alt = dt[is.na(aa_alt), which = T]
  if(length(invalid_alt) > 0) {
    problems[["invalid_aa_alt"]] = dt[invalid_alt, input_id]
    dt = dt[! invalid_alt]
    codons = codons[-invalid_alt]
  }
  
  dt[, alt_possible := mapply(function(codon, alt) length(unlist(codon_snvs_to_aa[[codon]][[alt]])) > 0, 
                                                 as.character(codons), aa_alt)]
  impossible_alt = dt[alt_possible == FALSE, which = T]
  if(length(impossible_alt) > 0) {
    problems[["aa_alt_not_possible"]] = dt[impossible_alt, input_id]
    dt = dt[! impossible_alt]
  }
  
  if(length(problems) == 0) {
    return(NULL)
  }
  return(problems)
}