#' Annotate variants
#' 
#' Annotates CESAnaysis MAF data with reference genome and gene data
#'
#' @param cesa CESAnalysis object
#' @export
annotate_variants <- function(cesa) {
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis", call. = F)
  }
  cesa = update_cesa_history(cesa, match.call())
  RefCDS = .ces_ref_data[[cesa@ref_key]]$RefCDS
  gr_genes = .ces_ref_data[[cesa@ref_key]]$gr_genes
  bsg = .ces_ref_data[[cesa@ref_key]]$genome
  cesa@advanced$annotated = TRUE
  
  # non-SNVs are not supported in selection functions yet, so not bothering to annotate them correctly
  # all non-SNV annotations will get set to NA at the end
  variants = unique(cesa@maf[, .(Chromosome, Start_Position, Reference_Allele, Tumor_Allele, variant_type)])
  
  # calling as.integer to avoid problem with scientific notation (100000->"1e05")
  variants[variant_type == "snv", snv_id := paste0(Chromosome, ':', as.integer(Start_Position), '_', Reference_Allele, '>', Tumor_Allele)]
  
  # Same-sample variants with 2bp of other variants get set aside as likely MNVs
  # MNVs are only possible in sample/chromosome combinations with more than one MAF record
  # Only new SNVs can be involved in MNVs
  # Will save the possible MNVs now to check later
  poss_mnv = cesa@maf[variant_type == 'snv']
  if ("variant_id" %in% names(cesa@maf)) {
    poss_mnv = poss_mnv[is.na(variant_id)] # restricts to previously unchecked variants
  }
  poss_mnv = poss_mnv[order(Unique_Patient_Identifier, 
                   Chromosome, Start_Position)][, .SD[.N > 1] , by = c("Unique_Patient_Identifier", "Chromosome")]
  
  # Don't need to re-annotate variants previously annotated in the CESAnalysis (on previous load_maf calls, usually)
  if(! is.null(cesa@mutations$snv)) {
    variants = variants[! cesa@mutations$snv$snv_id, on = "snv_id"]
    variants = variants[variant_type != 'likely_mnv'] # these have already been annotated, too
  }
  
  # There may be no variants left, especially if multiple CESAnalyses are being used on one set of variants
  # But, we can't just return, because MAF and coverage annotations must be updated
  variant_grs = GenomicRanges::GRanges()
  if(variants[, .N] > 0) {
    variant_grs = GenomicRanges::makeGRangesFromDataFrame(variants, seqnames.field = "Chromosome", 
                                                        start.field = "Start_Position", end.field = "Start_Position")
  }

  
  
  # Get gene annotations
  ## Note: If RefCDS/gr_genes are transcript-level, the ranges and annotations will refer to transcripts, not genes
  nearest = as.data.table(GenomicRanges::distanceToNearest(variant_grs, gr_genes, select = "all"))
  
  # convert the "subjectHits" index returned by the distanceToNearest function to the corresponding gene name
  gr_gene_names = GenomicRanges::mcols(gr_genes)["names"][,1]
  nearest[, gene := nearest[, gr_gene_names[subjectHits]]]
  
  # remove all but one of multiple hits from one record to the same gene
  # this happens when a record overlaps multiple exons in the reference data from different transcripts
  nearest = nearest[! duplicated(nearest[, .(queryHits, gene)])] 
  
  # queryHits column gives snv table row, gene gives associated gene
  # some records have multiple matching genes; combine them into variable-length vector within each snv table row
  # also grab the distance from the first matching gene (distances are always equal on multiple hits, since we asked for the nearest gene)
  genes_by_variant_row = nearest[, .(genes = list(gene), dist = distance[1]), by = "queryHits"]
  variants[, genes := genes_by_variant_row$genes]
  
  # handle edge case of just 1 SNV being annotated (which incorrectly makes genes non-list)
  if (is(variants$genes, "character")) {
    listed_genes = list(genes_by_variant_row$genes)
    variants[, genes := NULL]
    variants[, genes := listed_genes]
  }
   variants[, dist := genes_by_variant_row$dist]
   variants[, intergenic := dist > 0]
   
   snvs = variants[variant_type == "snv", -"variant_type"]
   indels = variants[variant_type == "indel"]

  
  # Intergenic records definitely aren't coding, so leave them out
  # some of these will still be non-coding and get filtered out later
  aac = data.table()
  if (snvs[dist == 0, .N] > 0) {
    aac = snvs[dist == 0, .(Chromosome, Start_Position, Reference_Allele, Tumor_Allele, gene = unlist(genes)), by = "snv_id"]
    genes_with_aac = unique(unlist(aac$gene))
  }

  # If RefCDS object has been created using this package's build_RefCDS function,
  # the presence of real_gene_name indicates that the object has one record per passing transcript rather that gene.
  # We'll handle this situation by making note of the true gene and transcript names.
  
  if (aac[, .N] > 0) {
    if (! is.null(RefCDS[[1]]$real_gene_name)) {
      # for efficiency, pull out all the necessary RefCDS data together, then split it up
      ref_subset = lapply(RefCDS[genes_with_aac], function(x) list(meta = list(strand = x$strand, cds = x$protein_id, gene = x$real_gene_name), 
                                                                   intervals_cds = x$intervals_cds))
      meta = rbindlist(lapply(lapply(ref_subset, function(x) x$meta), as.data.table))
      meta[, transcript_id := genes_with_aac]
      coding_ints = rbindlist(lapply(lapply(ref_subset, function(x) x$intervals_cds), as.data.table), idcol = "transcript_id")
      setnames(coding_ints, old = c("V1", "V2"), new = c("start", "end"))
      coding_ints = meta[coding_ints, on = "transcript_id" ]
    } else {
      ref_subset = lapply(RefCDS[genes_with_aac], function(x) list(meta = list(strand = x$strand, cds = x$protein_id, gene = x$gene_name), 
                                                                   intervals_cds = x$intervals_cds))
      coding_ints = rbindlist(lapply(lapply(ref_subset, function(x) x$intervals_cds), as.data.table), idcol = "gene")
      setnames(coding_ints, old = c("V1", "V2"), new = c("start", "end"))
      meta = rbindlist(lapply(lapply(ref_subset, function(x) x$meta), as.data.table))
      coding_ints = meta[coding_ints, on = "gene" ]
    } 
    # CDS intervals should be in genomic order within each cds
    coding_ints[strand == 1, cds_order :=  seq_len(.N), by = "cds"]
    coding_ints[strand == -1, cds_order := rev(seq_len(.N)), by = "cds"]
    coding_ints = coding_ints[order(cds, cds_order)]
    coding_ints[, cum_cds_width := cumsum(end - start + 1), by = "cds" ]
    
    aac[, tmp_end_pos := Start_Position] # need two columns for foverlaps
    setkey(aac, "gene", "Start_Position", "tmp_end_pos") # yeah, it's called gene for transcript_id currently, too
    
    ## To-do: switch gene/transcript_id here
    cds_hits = foverlaps(coding_ints, aac, by.x = c("gene", "start", "end"), type = "any", nomatch = NULL)
  } else {
    cds_hits = data.table()
  }

  if (cds_hits[, .N] > 0) {
    cds_hits[strand == -1, nt_pos := cum_cds_width - (Start_Position - start)]
    cds_hits[strand == 1, nt_pos := cum_cds_width - (end - Start_Position)]
    cds_hits[, aa_pos := ceiling(nt_pos / 3)]
    cds_hits[, codon_pos := nt_pos %% 3 ] # note this gives 0-based codon index 
    cds_hits[codon_pos == 0, codon_pos := 3]
    
    cds_genes = unique(cds_hits$gene)
    seqs = DNAStringSet(lapply(cds_genes, function(x) RefCDS[[x]]$seq_cds))
    names(seqs) = cds_genes
    ref_seqs = Biostrings::subseq(seqs[cds_hits$gene], start = cds_hits$nt_pos - cds_hits$codon_pos + 1, width = 3)
    
    # Records that don't overlap a CDS are not coding mutations
    aac = cds_hits
    
    # For efficiency, convert to character, and then substitute in the tumor allele and get aa_alt
    aa_ref = as.character(Biostrings::translate(ref_seqs, no.init.codon = T), use.names = F)
    aac_alt_allele = Biostrings::DNAStringSet(cds_hits$Tumor_Allele)
    aac_alt_allele[cds_hits$strand == -1] = Biostrings::complement(aac_alt_allele[cds_hits$strand == -1])
    alt_seqs = as.character(ref_seqs)
    substr(alt_seqs, cds_hits$codon_pos, cds_hits$codon_pos) = as.character(aac_alt_allele)
    aa_alt = as.character(Biostrings::translate(DNAStringSet(alt_seqs), no.init.codon = T))
    
    aac[, aachange := paste0(aa_ref, aa_pos, aa_alt)]
    aac[, aa_ref := seqinr::aaa(aa_ref)]
    aac[aa_ref == "Stp", aa_ref := "STOP"]
    aac[, aa_alt := seqinr::aaa(aa_alt)]
    aac[aa_alt == "Stp", aa_alt := "STOP"]
    
    setnames(aac, old = c("cds", "Chromosome", "Start_Position"), new = c("pid", "chr", "pos"))
    aac[, aac_id := paste0(gene, "_", aachange, "_", pid)]
    aac[, c("start", "end", "cds_order", "cum_cds_width") := NULL]
    
    # Some AACs will come from >1 distinct SNV in the MAF data; only need one of each for annotation purposes
    aac = unique(aac, by = "aac_id")
    
    # To make things easier later, determine which AAC records are possibly spanning splice sites
    # Genomically, all bases must be within 3bp of a splice position
    # Later, we'll do a more careful annotation of splice site distance
    # For now, get CDS intervals for relevant genes; potential splice sites will be included in the start/end positions of these
    approx_splice_sites = unique(rbindlist(lapply(RefCDS[genes_with_aac], function(x) return(list(chr = x$chr, splice_pos = x$intervals_cds)))))
    approx_splice_sites[, start := splice_pos - 3]
    approx_splice_sites[, end := splice_pos + 3]
    setkey(aac, "chr", "pos", "tmp_end_pos")
    possible_splice_aac_id = foverlaps(approx_splice_sites, aac, by.x = c("chr", "start", "end"), type = "any", nomatch = NULL)[, aac_id]
    setkey(aac, "aac_id")
    aac[, possible_splice := F][possible_splice_aac_id, possible_splice := T]
    
    
    # function to get reference positions of every nucleotide in a codon (needed for near splice sites)
    calc_ref_pos = function(gene, nt_pos, strand) {
      intervals = RefCDS[[gene]]$intervals_cds
      pos = apply(intervals, 1, function(x) abs(x[1]-x[2]) + 1)
      if(strand == -1) {
        pos = rev(pos)
        if (nrow(intervals) > 1) {
          intervals = apply(intervals, 2, rev) # because 1-row matrix would get coerced to vector by apply
        }
      }
      
      nt_so_far = cumsum(pos)
      ind = match(TRUE, nt_so_far >= nt_pos)
      
      previous_bases = ifelse(ind == 1, 0, nt_so_far[ind - 1])
      if (strand == 1) {
        ref_pos = intervals[ind, 1] + nt_pos - previous_bases - 1
      } else {
        if(ncol(intervals) < 2 || nrow(intervals) < ind) {
          return(0)
        }
        ref_pos = intervals[ind, 2] - nt_pos + previous_bases + 1
        
      }
      return(ref_pos)
    }
    
    # away from splice sites, all codon nt positions are adjacent
    aac[possible_splice == F, nt1_pos := pos - strand * (codon_pos - 1)]
    aac[possible_splice == F, nt2_pos := pos - strand * (codon_pos - 2)]
    aac[possible_splice == F, nt3_pos := pos - strand * (codon_pos - 3)]
    
    # figure out nt positions for codons that might span splice sites
    aac[possible_splice == T, nt1_pos := mapply(calc_ref_pos, gene, nt_pos - codon_pos + 1, strand)]
    aac[possible_splice == T, nt2_pos := mapply(calc_ref_pos, gene, nt_pos - codon_pos + 2, strand)]
    aac[possible_splice == T, nt3_pos := mapply(calc_ref_pos, gene, nt_pos - codon_pos + 3, strand)]
    
    # get 5'->3' (i.e., standard genomic order) trinucleotide contexts
    first_triseq = BSgenome::getSeq(bsg, names = aac$chr, start =  aac$nt1_pos - 1, end =  aac$nt1_pos + 1)
    second_triseq = BSgenome::getSeq(bsg, names = aac$chr, start =  aac$nt2_pos - 1, end =  aac$nt2_pos + 1)
    third_triseq = BSgenome::getSeq(bsg, names = aac$chr, start =  aac$nt3_pos - 1, end =  aac$nt3_pos + 1)
    
    # put together the three reference nucleotides of each codon (in coding order, not genomic order)
    coding_sequences = Biostrings::xscat(Biostrings::subseq(first_triseq, start = 2, width = 1),
                                         Biostrings::subseq(second_triseq, start = 2, width = 1),
                                         Biostrings::subseq(third_triseq, start = 2, width = 1))
    
    aac[, nt1_ref := as.character(Biostrings::subseq(coding_sequences, start = 1, width = 1))]
    aac[, nt2_ref := as.character(Biostrings::subseq(coding_sequences, start = 2, width = 1))]
    aac[, nt3_ref := as.character(Biostrings::subseq(coding_sequences, start = 3, width = 1))]
    
    
    # take complement of negative-strand coding sequences (but no need to reverse)
    coding_sequences[aac$strand == -1] = Biostrings::complement(coding_sequences[aac$strand == -1])
    aac[, coding_seq := as.character(coding_sequences)]
    
    # for each aa_mut, get a three-item list representing SNVs in first, second, third positions of the reference codon
    # that cause the same amino acid change as aa_alt (which will always include the ntchange that came from the MAF)
    possible_snvs = mapply(function(x, y) codon_snvs_to_aa[[x]][[y]], as.character(coding_sequences), aac$aa_alt, SIMPLIFY = F)
    
    
    nt1_snvs = lapply(possible_snvs, function(x) x[[1]])
    aac[, nt1_snvs := nt1_snvs]
    
    nt2_snvs = lapply(possible_snvs, function(x) x[[2]])
    aac[, nt2_snvs := nt2_snvs]
    
    nt3_snvs = lapply(possible_snvs, function(x) x[[3]])
    aac[, nt3_snvs := nt3_snvs]
    
    # create table of all snv_ids, then fix strand of alt alleles
    nt1_dt = aac[, .(chr, strand, pos = nt1_pos, ref = nt1_ref, alt = unlist(nt1_snvs)), by = "aac_id"]
    nt2_dt = aac[, .(chr, strand, pos = nt2_pos, ref = nt2_ref, alt = unlist(nt2_snvs)), by = "aac_id"]
    nt3_dt = aac[, .(chr, strand, pos = nt3_pos, ref = nt3_ref, alt = unlist(nt3_snvs)), by = "aac_id"]
    snv_table = rbind(nt1_dt, nt2_dt, nt3_dt)
    snv_table = snv_table[! is.na(alt)] # drop positions where there is no possible SNV that causes the amino acid change of interest
    snv_table[strand == -1, alt := seqinr::comp(alt, forceToLower = F)]
    snv_table[, strand := NULL] # SNVs don't have strandedness
    snv_table[, snv_id := paste0(chr, ':', pos, "_", ref, '>', alt)]
    # Merge in the gene annotations from the original snvs table
    # There will be no gene annotations for the SNVs that are part of AACs but
    # not present in the MAF data, so we'll have to do annotation again on those
    snv_table = snvs[, .(snv_id, genes, intergenic)][snv_table, , on = "snv_id"]
    
    # annotate amino acid mutations with all associated SNVs
    snvs_by_aa_mut = snv_table[, .(constituent_snvs = list(snv_id)), by = "aac_id"]
    aac[snvs_by_aa_mut, constituent_snvs := constituent_snvs, on = "aac_id"]
    snv_table = snv_table[, .(chr, pos, ref, alt, genes, intergenic, assoc_aac = list(sort(aac_id))), by = "snv_id"]
    
    # add noncoding SNVs to SNV table
    noncoding = snvs[! snv_id %in% unlist(aac$constituent_snvs)]
    noncoding = noncoding[, .(chr = Chromosome, pos = Start_Position, ref = Reference_Allele, alt = Tumor_Allele,
                              genes, intergenic)]
    noncoding[, snv_id := paste0(chr, ':', pos, "_", ref, '>', alt)]
    setcolorder(noncoding, "snv_id")
    # dt can't handle list columns in unique call, but if it could result would be the same
    noncoding = unique(noncoding, by = c("chr", "pos", "ref", "alt")) 
    noncoding[, assoc_aac := list(as.list(NA_character_))] # to match coding records
    snv_table = rbind(snv_table, noncoding)
    
  } else {
    # When there are no CDS hits, all the SNVs are noncoding
    snv_table = snvs[, .(chr = Chromosome, pos = Start_Position, ref = Reference_Allele, alt = Tumor_Allele,
                         genes = genes, intergenic, assoc_aac = list(NA_character_)), by = "snv_id"]
    aac = aac[0] # empty the table
  }
  
	# uniquify
  snv_table = unique(snv_table, by = "snv_id")
  
  # Get gene annotations for those that don't have it yet (identified by intergenic == NA)
  # See comments from round 1
  snvs_needing_anno = snv_table[is.na(intergenic)]
  if (snvs_needing_anno[, .N] > 0) {
    snv_gr_round2 = GenomicRanges::makeGRangesFromDataFrame(snvs_needing_anno,  seqnames.field = "chr", 
                                                            start.field = "pos", end.field = "pos")
    
    nearest = as.data.table(GenomicRanges::distanceToNearest(snv_gr_round2, gr_genes, select = "all"))
    possible_genes = GenomicRanges::mcols(gr_genes)["names"][,1]
    nearest[, gene := nearest[, possible_genes[subjectHits]]]
    nearest = nearest[! duplicated(nearest[, .(queryHits, gene)])] 
    genes_by_variant_row = nearest[, .(genes = list(gene), dist = distance[1]), by = "queryHits"]
    snvs_needing_anno[, genes := genes_by_variant_row$genes]
    snvs_needing_anno[, intergenic := genes_by_variant_row$dist > 0]
    
    setkey(snv_table, "snv_id")
    snv_table[snvs_needing_anno$snv_id, intergenic := snvs_needing_anno$intergenic]
    snv_table[snvs_needing_anno$snv_id, genes := list(snvs_needing_anno$genes)]
  }
  

  #  get deconstructSigs-style trinuc context of each SNV ID
  if(snv_table[, .N] > 0) {
    genomic_context = BSgenome::getSeq(bsg, snv_table$chr,start=snv_table$pos - 1,
                                       end=snv_table$pos + 1,
                                       as.character = TRUE)
    trinuc_mut_ids = paste0(genomic_context,":", snv_table$alt)
    
    # deconstructSigs_notations is a keyed table in CES sysdata
    snv_table[, trinuc_mut := deconstructSigs_notations[.(genomic_context, snv_table$alt), deconstructSigs_ID]] 
  }

  
  maf = cesa@maf
  maf[variant_type == "snv", variant_id := paste0(Chromosome, ':', Start_Position, '_', Reference_Allele, '>', Tumor_Allele)]
  
  # clean up aac table, except when it's empty
  if (aac[, .N] > 0) {
    aac_table = aac
    aac_table = aac_table[, .(aac_id, chr, gene, strand, pid, aachange, aa_ref, aa_pos, aa_alt, possible_splice, nt1_pos, nt2_pos, nt3_pos, coding_seq, constituent_snvs)]
    setcolorder(aac_table, c("aac_id", "gene", "aachange", "strand"))
    
    # If any trinucleotide mutation comes up NA--usually due to an ambiguous N in the genomic trinuc context--remove record from analysis
    bad_trinuc_context = which(is.na(snv_table$trinuc_mut))
    num_bad = length(bad_trinuc_context)
    if (num_bad > 0) {
      bad_ids = snv_table[bad_trinuc_context, snv_id]
      bad_trinuc_context_maf <- maf[variant_id %in% bad_ids, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
      maf <- maf[! variant_id %in% bad_ids]
      message(paste("Note:", num_bad, "MAF records were excluded from analysis because of ambiguous trinucleotide context (likely N's in the reference genome)."))
      bad_trinuc_context_maf$Exclusion_Reason = "ambiguous_trinuc_context"
      cesa@excluded = rbind(cesa@excluded, bad_trinuc_context_maf)
      
      # remove the bad record from SNV and AAC tables
      bad_aa = unlist(snv_table[bad_trinuc_context, assoc_aac])
      snv_table = snv_table[! bad_trinuc_context]
      aac_table = aac_table[! aac_id %in% bad_aa]
    }
  }
  
	# Annotate SNVs that are at essential splice sites accroding to RefCDS, and then apply to any associated AACs
	# Note that we're not keeping track of which genes these splice site positions apply to
	genes_in_data = unique(unlist(snv_table$genes))
	ind_no_splice = sapply(RefCDS[genes_in_data], function(x) length(x$intervals_splice) == 0)
	genes_with_splice_sites = genes_in_data[! ind_no_splice]
	snv_table[, essential_splice := F]
	if (length(genes_with_splice_sites) > 0) {
	  essential_splice_sites = unique(rbindlist(lapply(RefCDS[genes_with_splice_sites], function(x) return(list(chr = x$chr, start = x$intervals_splice)))))
	  essential_splice_sites[, end := start]
	  snv_table[, tmp_end_pos := pos]
	  setkey(snv_table, "chr", "pos", "tmp_end_pos")
	  essential_splice_snv_id = foverlaps(essential_splice_sites, snv_table, by.x = c("chr", "start", "end"), type = "any", nomatch = NULL)[, snv_id]
	  snv_table[, tmp_end_pos := NULL]
	  setkey(snv_table, "snv_id")
	  snv_table[essential_splice_snv_id, essential_splice := T]
	}

	
	# If an AAC has any constituent SNV at an essential splice site, we'll say the AAC is at an essential splice site
	if (aac[, .N] > 0) {
	  setkey(aac_table, "aac_id")
	  essential_splice_aac_id = snv_table[essential_splice == T & ! is.na(assoc_aac), unique(unlist(assoc_aac))]
	  aac_table[, essential_splice := F][essential_splice_aac_id, essential_splice := T]
	  aac_table[, possible_splice := NULL]
	}

	
	if (aac[, .N] > 0) {
	  cesa@mutations[["amino_acid_change"]] = unique(rbind(cesa@mutations$amino_acid_change, aac_table, fill = T), by = "aac_id")
	  setkey(cesa@mutations$amino_acid_change, "aac_id")
	}
  cesa@mutations[["snv"]] = unique(rbind(cesa@mutations$snv, snv_table, fill = T), by = "snv_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  # clear fields that come from annotate_variants (if previously any annotations)
  if (! is.null(maf$genes)) {
    maf[, c("genes", "assoc_aac") := NULL]
  }
  
  # add genes and assoc_aac to MAF records
  # use of _tmp names required as of data.table 1.13.2 to keep join from failing
  maf[, c("genes_tmp", "assoc_aac_tmp") := list(list(NA_character_), list(NA_character_))]
  maf[cesa@mutations$snv, c("genes_tmp", "assoc_aac_tmp") := list(genes, assoc_aac), on = c(variant_id = "snv_id")]
  intergenic_snv = cesa@mutations$snv[intergenic == T, snv_id]
  genic_indels = indels[dist == 0]
  maf[genic_indels, c("genes_tmp", "assoc_aac_tmp") := list(genes, list(NA_character_)), 
      on = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Allele")]
  maf[intergenic_snv, genes_tmp := list(NA_character_), on= "variant_id"]
  setnames(maf, c("genes_tmp", "assoc_aac_tmp"), c("genes", "assoc_aac"))
  setcolorder(maf, colnames(cesa@maf)) # put original columns back in the front
  
  # Handle the possible MNVs identified early in the function
  if (poss_mnv[, .N] > 0) {
    poss_mnv[, dist_to_prev := c(Inf, diff(Start_Position)), by = c("Unique_Patient_Identifier", "Chromosome")]
    poss_mnv[, dist_to_next := c(dist_to_prev[2:.N], Inf), by = c("Unique_Patient_Identifier", "Chromosome")]
    poss_mnv[dist_to_prev < 3 | dist_to_next < 3, is_mnv := T]
    maf[poss_mnv, is_mnv := is_mnv, on = c("Unique_Patient_Identifier", "Chromosome", "Start_Position")]
    num_mnv = maf[is_mnv == T, .N]
    
    # When called internally, skip message (probably called by add_variants)
    if (num_mnv > 0 & cesa@advanced$recording == TRUE) {
      msg = paste0("Note: ", num_mnv, " MAF records ",
                   "are within 2 bp of other same-sample variants and will be classified as \"likely_mnv\" instead of snv/indel. ",
                   "(That is, they're probably multi-nucleotide variants arising from single events.)")
      pretty_message(msg)
    }
    maf[is_mnv == TRUE, c("variant_type", "variant_id", "assoc_aac") := list("likely_mnv", NA_character_, list(NA_character_))] 
    maf[, is_mnv := NULL]
  }
  cesa@maf = maf
  
  # Add covered_in to annotation tables and return
  # Internal note: Confusingly, update_covered_in also has a side effect of updating
  # cached output of select_variants, because whether using annotate_variants,
  # add_covered_regions, or add_variants, it always gets called at the end.
  return(update_covered_in(cesa))
	
}

