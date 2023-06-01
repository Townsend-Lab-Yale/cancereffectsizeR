#' Annotate SNVs
#' 
#' Annotates SNVs with reference information, including protein changes.
#'
#' @param refset CES reference data set (e.g., from the ces.refset.hg38 data package).
#' @param variants MAF-like data.table of SNVs.
#' @keywords internal
annotate_snvs = function(snvs, refset) {
  
  output = list(amino_acid_change = aac_annotation_template, snv = snv_annotation_template,
                aac_snv_key = aac_snv_key_template)
  if(snvs[, .N] == 0) {
    return(output)
  }
  RefCDS = refset$RefCDS
  gr_cds = refset$gr_genes
  bsg = refset$genome
  aac_snv_key = copy(aac_snv_key_template)
  
  snvs = get_nearest_refset_entries(variants = snvs, refset = refset)
  setnames(snvs, 'variant_id', 'snv_id')
  
  # Intergenic records definitely aren't coding, so leave them out
  # some of these will still be non-coding and get filtered out later
  aac = data.table()
  if (snvs[dist == 0, .N] > 0) {
    aac = snvs[dist == 0, .(Chromosome, Start_Position, Reference_Allele, Tumor_Allele, refcds_entry = unlist(cds)), by = "snv_id"]
    entries_with_aac = unique(unlist(aac$refcds_entry))
  }
  
  # If RefCDS object has been created using this package's build_RefCDS function,
  # the presence of real_gene_name indicates that the object has one record per passing transcript rather that gene.
  # We'll handle this situation by making note of the true gene and transcript names.
  cds_hits = data.table()
  if (aac[, .N] > 0) {
    # for efficiency, pull out all the necessary RefCDS data together, then split it up
    if(is.null(RefCDS[[1]]$real_gene_name)) {
      # real_gene_name is used for pid-based RefCDS objects. When null, gene name is already correct
      ref_subset = lapply(RefCDS[entries_with_aac], function(x) list(meta = list(strand = x$strand, cds = x$protein_id, entry_name = x$gene_name,
                                                                                 gene_name = x$gene_name), 
                                                                     intervals_cds = x$intervals_cds))
    } else {
      ref_subset = lapply(RefCDS[entries_with_aac], function(x) list(meta = list(strand = x$strand, cds = x$protein_id, entry_name = x$gene_name,
                                                                                 gene_name = x$real_gene_name), 
                                                                     intervals_cds = x$intervals_cds))
    }
    meta = rbindlist(lapply(ref_subset, function(x) x$meta))
    coding_ints = rbindlist(lapply(lapply(ref_subset, function(x) x$intervals_cds), as.data.table), idcol = "entry_name")
    setnames(coding_ints, old = c("V1", "V2"), new = c("start", "end"))
    coding_ints = meta[coding_ints, on = "entry_name"]
    
    # CDS intervals should be in genomic order within each cds
    coding_ints[strand == 1, cds_order := seq_len(.N), by = "cds"]
    coding_ints[strand == -1, cds_order := rev(seq_len(.N)), by = "cds"]
    coding_ints = coding_ints[order(cds, cds_order)]
    coding_ints[, cum_cds_width := cumsum(end - start + 1), by = "cds"]
    
    aac[, tmp_end_pos := Start_Position] # need two columns for foverlaps (as of now)
    setkey(aac, "refcds_entry", "Start_Position", "tmp_end_pos") # yeah, it's called gene for transcript_id currently, too
    
    cds_hits = foverlaps(coding_ints, aac, by.x = c("entry_name", "start", "end"), type = "any", nomatch = NULL)
  }
  
  if (cds_hits[, .N] > 0) {
    cds_hits[strand == -1, nt_pos := cum_cds_width - (Start_Position - start)]
    cds_hits[strand == 1, nt_pos := cum_cds_width - (end - Start_Position)]
    cds_hits[, aa_pos := ceiling(nt_pos / 3)]
    cds_hits[, codon_pos := nt_pos %% 3 ] # note this gives 0-based codon index 
    cds_hits[codon_pos == 0, codon_pos := 3]
    
    cds_names = unique(cds_hits$entry_name)
    seqs = DNAStringSet(lapply(RefCDS[cds_names], '[[', 'seq_cds'))
    ref_seqs = Biostrings::subseq(seqs[cds_hits$entry_name], start = cds_hits$nt_pos - cds_hits$codon_pos + 1, width = 3)
    
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
    aac[, aac_id := paste0(gene_name, "_", aachange, "_", pid)]
    aac[, c("start", "end", "cds_order", "cum_cds_width") := NULL]
    
    # Some AACs will come from >1 distinct SNV in the MAF data; only need one of each for annotation purposes
    aac = unique(aac, by = "aac_id")
    
    # To make things easier later, determine which AAC records are possibly spanning splice sites
    # Genomically, all bases must be within 2bp of a splice position (or codon finishes before crossing site)
    # Later, we'll do a more careful annotation of splice site distance
    # For now, get CDS intervals for relevant genes; potential splice sites will be included in the start/end positions of these
    approx_splice_sites = unique(rbindlist(lapply(RefCDS[entries_with_aac], function(x) return(list(chr = x$chr, splice_pos = x$intervals_cds)))))
    approx_splice_sites[, start := splice_pos - 2]
    approx_splice_sites[, end := splice_pos + 2]
    setkey(aac, "chr", "pos", "tmp_end_pos")
    possible_splice_aac_id = foverlaps(approx_splice_sites, aac, by.x = c("chr", "start", "end"), type = "any", nomatch = NULL)[, aac_id]
    setkey(aac, "aac_id")
    aac[, possible_splice := F][possible_splice_aac_id, possible_splice := T]
    
    # figure out nt positions for codons that might span splice sites
    poss_splice_ints = coding_ints[aac[possible_splice == TRUE, unique(entry_name)], on = 'entry_name']
    single_exon_entries = poss_splice_ints[, .N, by = "entry_name"][N == 1, entry_name]
    poss_splice_ints = poss_splice_ints[! single_exon_entries, on = 'entry_name']
    aac[single_exon_entries, possible_splice := FALSE, on = 'entry_name']
    
    # away from splice sites, all codon nt positions are adjacent
    aac[possible_splice == F, nt1_pos := pos - strand * (codon_pos - 1)]
    aac[possible_splice == F, nt2_pos := pos - strand * (codon_pos - 2)]
    aac[possible_splice == F, nt3_pos := pos - strand * (codon_pos - 3)]
    
    
    if (poss_splice_ints[, .N] > 0) {
      poss_splice_ints[, exon_length := abs(end - start) + 1]
      poss_splice_ints[, previous_bases := c(0, cumsum(exon_length[1:(.N-1)])), by = 'entry_name']
      poss_splice_ints[strand == 1, ref_offset := start - previous_bases - 1]
      poss_splice_ints[strand == -1, ref_offset := end + previous_bases + 1]
      
      # some SNVs may have multiple entries; splice status may vary
      poss_splice_table = setDT(aac[possible_splice == T, .(tmp_nt_pos = nt_pos - codon_pos + 1:3, tmp_codon_pos = 1:3, strand), by = c('snv_id', 'entry_name')])
      ref_offsets = poss_splice_ints[poss_splice_table, ref_offset, on = .(entry_name, cum_cds_width >= tmp_nt_pos), by = .EACHI, mult = "first"]
      poss_splice_table[, nt_ref_pos := ref_offsets$ref_offset + tmp_nt_pos * strand]
      aac[poss_splice_table[tmp_codon_pos == 1], nt1_pos := nt_ref_pos, on = c('snv_id', 'entry_name')]
      aac[poss_splice_table[tmp_codon_pos == 2], nt2_pos := nt_ref_pos, on = c('snv_id', 'entry_name')]
      aac[poss_splice_table[tmp_codon_pos == 3], nt3_pos := nt_ref_pos, on = c('snv_id', 'entry_name')]
    }
    
    
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
    aac$nt1_snvs = nt1_snvs
    
    nt2_snvs = lapply(possible_snvs, function(x) x[[2]])
    aac$nt2_snvs = nt2_snvs
    
    nt3_snvs = lapply(possible_snvs, function(x) x[[3]])
    aac$nt3_snvs = nt3_snvs
    
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
    new_snvs = setdiff(snv_table$snv_id, snvs$snv_id)
    if(length(new_snvs) > 0) {
      new_snv_anno = get_nearest_refset_entries(variants = snv_table[new_snvs,
                                                                     .(Chromosome = chr, Start_Position = pos,
                                                                       Reference_Allele = ref, Tumor_Allele = alt, variant_id = snv_id,
                                                                       variant_type = 'snv'), on = 'snv_id'], 
                                                refset = refset)
      setnames(new_snv_anno, 'variant_id', 'snv_id')
      snvs = rbind(snvs, new_snv_anno)
    }
    snv_table[snvs, c("genes", "nearest_pid", "cds", "dist") := .(genes, nearest_pid, cds, dist), on = 'snv_id']
    
    # annotate amino acid mutations with all associated SNVs
    snvs_by_aa_mut = snv_table[, .(constituent_snvs = list(snv_id)), by = "aac_id"]
    aac[snvs_by_aa_mut, constituent_snvs := constituent_snvs, on = "aac_id"]
    
    aac_snv_key = snv_table[, .(aac_id, snv_id)]
    aac_snv_key[, multi_anno_site := uniqueN(aac_id) > 1, by = 'snv_id']
    snv_table = snv_table[, .(snv_id, chr, pos, ref, alt, genes, nearest_pid, cds, dist)] # gets uniquified shortly
    
    # add noncoding SNVs to SNV table
    noncoding = snvs[! unlist(aac$constituent_snvs), on = 'snv_id']
    noncoding = setDT(noncoding[, .(chr = Chromosome, pos = Start_Position, ref = Reference_Allele, alt = Tumor_Allele,
                                    genes, nearest_pid, dist, cds)])
    noncoding[, snv_id := paste0(chr, ':', pos, "_", ref, '>', alt)]
    setcolorder(noncoding, "snv_id")
    # dt can't handle list columns in unique call, but if it could result would be the same
    noncoding = unique(noncoding, by = c("chr", "pos", "ref", "alt")) 
    snv_table = rbind(snv_table, noncoding)
    
  } else {
    # When there are no CDS hits, all the SNVs are noncoding
    snv_table = snvs[, .(chr = Chromosome, pos = Start_Position, ref = Reference_Allele, alt = Tumor_Allele,
                         genes, nearest_pid, cds, dist), by = "snv_id"]
    aac = aac[0] # empty the table
    aac[, aac_id := character()]
  }
  
  # Remove redundant rows (can get created when multiple AACs have shared constituent SNVs).
  snv_table = unique(snv_table, by = "snv_id")
  snv_table[, intergenic := dist > 0] # Not really intergenic, just not coding (need to fix this)
  
  genomic_context = BSgenome::getSeq(bsg, snv_table$chr,start=snv_table$pos - 1,
                                     end=snv_table$pos + 1,
                                     as.character = TRUE)
  trinuc_mut_ids = paste0(genomic_context,":", snv_table$alt)
  
  # deconstructSigs_notations is a keyed table in CES sysdata
  snv_table[, trinuc_mut := deconstructSigs_notations[.(genomic_context, snv_table$alt), deconstructSigs_ID]]
  
  # Annotate SNVs that are at essential splice sites according to RefCDS, and then apply to any associated AACs.
  # Note that we're not keeping track of which genes these splice site positions apply to.
  cds_in_data = unique(unlist(snv_table$cds))
  ind_no_splice = sapply(RefCDS[cds_in_data], function(x) length(x$intervals_splice) == 0)
  cds_with_splice_sites = cds_in_data[! ind_no_splice]
  snv_table[, essential_splice := F]
  aac[, essential_splice := F]
  if (length(cds_with_splice_sites) > 0) {
    essential_splice_sites = unique(rbindlist(lapply(RefCDS[cds_with_splice_sites], function(x) return(list(chr = x$chr, start = x$intervals_splice)))))
    essential_splice_sites[, end := start]
    snv_table[, tmp_end_pos := pos]
    setkey(snv_table, "chr", "pos", "tmp_end_pos")
    essential_splice_snv_id = foverlaps(essential_splice_sites, snv_table, by.x = c("chr", "start", "end"), type = "any", nomatch = NULL)[, snv_id]
    snv_table[, tmp_end_pos := NULL]
    setkey(snv_table, "snv_id")
    snv_table[essential_splice_snv_id, essential_splice := T]
    
    # Most essential splice sites are noncoding, but a few overlap coding regions
    essential_splice_aac_id = aac_snv_key[essential_splice_snv_id, unique(aac_id), on = 'snv_id', nomatch = NULL]
    aac[essential_splice_aac_id, essential_splice := T, on = 'aac_id']
  }
  
  # clean up aac table, except when it's empty
  if (aac[, .N] > 0) {
    # to do: eventually, probably want to keep the entry_name field (call it refcds_entry?)
    aac_table = aac[, .(aac_id, chr, gene = gene_name, strand, pid, aachange, aa_ref, aa_pos, aa_alt, nt1_pos, nt2_pos, nt3_pos, 
                        coding_seq, constituent_snvs, essential_splice)]
    setkey(aac_table, 'aac_id')
    setcolorder(aac_table, c("aac_id", "gene", "aachange", "strand"))
  } else {
    aac_table = data.table()
  }
  snv_table[, c("dist", "cds") := NULL]
  setkey(snv_table, "snv_id")
  
  return(list(amino_acid_change = aac_table, snv = snv_table, aac_snv_key = aac_snv_key))
}


