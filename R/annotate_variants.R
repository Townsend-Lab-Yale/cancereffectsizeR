#' Annotate variants with data from genome and reference transcripts
#' 
#' Annotates CESAnaysis MAF data with reference genome and gene data, keeping assignments consistent with dNdScv when possible
#'
#' @importFrom IRanges "%within%"
#' @param cesa CESAnalysis object
#' @export
annotate_variants <- function(cesa) {
  
  RefCDS = get_genome_data(cesa, "RefCDS")
  gr_genes = get_genome_data(cesa, "gr_genes")
  
  
  # non-SNVs are not supported in selection functions yet, so not bothering to annotate them correctly
  # all non-SNV annotations will get set to NA at the end of this function
  MAF = cesa@maf
	dndscv_gene_names = cesa@mutrates$gene
	dndscv_out_list = cesa@dndscv_out_list

	# Select just the reference genes that are in the data output from dNdScv
	is_in_dndscv = (GenomicRanges::mcols(gr_genes)["names"][,1] %in% dndscv_gene_names)
	gr_genes_in_data = gr_genes[is_in_dndscv]
	
	# assign SNV IDs to all SNVs
	MAF[, snv_id := paste0(Chromosome, ':', Start_Position, '_', Reference_Allele, '>', Tumor_Allele)]

	# Get mutation annotations from dNdScv and subset to SNVs
	# Note: Mutations that are tossed by dNdScv--such as those exclusive to the most hypermutated samples or likely ketaegis--
	#       probably don't appear in annotations. This is okay for now, but someday this may be addressed.
	dndscvout_annotref <- rbindlist(lapply(dndscv_out_list, function(x) x$annotmuts))
	dndscvout_annotref <- dndscvout_annotref[ref %in% c("A", "C", "G", "T") & mut %in% c("A", "C", "G", "T")]
	dndscv_coding_anno = unique(dndscvout_annotref[impact != "Essential_Splice", .(chr, pos, ref, mut, gene, pid, strand, aachange, ntchange, impact)])
	setnames(dndscv_coding_anno, "mut", "alt")
	
	# build table of amino acid mutations by starting with SNVs with gene annotation
	aac_table = unique(MAF[, .(chr = Chromosome, pos = Start_Position, ref = Reference_Allele, alt = Tumor_Allele), by = "snv_id"])
	
	# merge in coding annotations (dropping non-coding records)
	aac_table = merge.data.table(aac_table, dndscv_coding_anno, by = c("chr", "pos", "ref", "alt"), sort = F)
	aac_table[, aac_id := paste0(gene, "_", aachange, "_", pid)]
	aac_table = unique(aac_table, by = "aac_id") # will figure out associated snv_ids later
	
	aac_table[, aa_ref := seqinr::aaa(stringr::str_sub(aachange, start = 1, end = 1))]
	aac_table[aa_ref == "Stp", aa_ref := "STOP"] # seqinr returns Stp for stop codons; we want STOP
	aac_table[, aa_alt := seqinr::aaa(stringr::str_sub(aachange, start = -1, end = -1))]
	aac_table[aa_alt == "Stp", aa_alt := "STOP"]
	aac_table[, aa_pos := as.numeric(stringr::str_sub(aachange, start = 2, end = -2))]
	
	
	# to make things easier later, determine now which MAF records are near splice sites (within 3 bp)
	# First, get DS intervals for every gene; potential splice sites are all the start/end positions of these
	splice_sites = rbindlist(lapply(RefCDS, function(x) return(list(gene = x$gene_name, splice_pos = x$intervals_cds))))
	comb = merge.data.table(aac_table[,.(aac_id, gene, pos)], splice_sites, allow.cartesian = T)
	comb[, diff:= abs(splice_pos - pos)]
	comb = comb[, .(next_to_splice = any(diff <= 3)), by = aac_id]
	aac_table[comb, next_to_splice := next_to_splice, on = "aac_id"]
	
	# for debug, ensure dndscv's reference nt matches the MAF ref column
	# note that dNdScv gives the nucleotide on the gene's strand (which is odd, I think)
	aac_table[, nt_ref := stringr::str_sub(ntchange, start = 1, end = 1)]
  aac_table[strand == -1, nt_ref := seqinr::comp(nt_ref, forceToLower = F)]
	if (! aac_table[, identical(nt_ref, ref)]) {
	  stop(paste0("dNdScv's reference allele (from RefCDS) doesn't match MAF reference allele.\n",
	              "If you're using custom transcript data, it may not match the reference genome.\n",
	              "If you're using cancereffectsizeR reference data, please submit a bug report."), call. = F)
	}
	aac_table[, cds_nt_pos := as.numeric(stringr::str_sub(ntchange, start = 2, end = -2))]
	aac_table[, codon_pos := cds_nt_pos %% 3]
	aac_table[codon_pos==0, codon_pos := 3] # "codon 0" obtained from  (nuc_position % 3) is actually nt 3 of codon
	
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
	aac_table[next_to_splice == F, nt1_pos := pos - strand * (codon_pos - 1)]
	aac_table[next_to_splice == F, nt2_pos := pos - strand * (codon_pos - 2)]
	aac_table[next_to_splice == F, nt3_pos := pos - strand * (codon_pos - 3)]
	
	# figure out nt positions for codons that might span splice sites
	aac_table[next_to_splice == T, nt1_pos := mapply(calc_ref_pos, gene, cds_nt_pos - codon_pos + 1, strand)]
	aac_table[next_to_splice == T, nt2_pos := mapply(calc_ref_pos, gene, cds_nt_pos - codon_pos + 2, strand)]
	aac_table[next_to_splice == T, nt3_pos := mapply(calc_ref_pos, gene, cds_nt_pos - codon_pos + 3, strand)]
	
	# get 5'->3' (i.e., standard genomic order) trinucleotide contexts
	first_triseq = BSgenome::getSeq(cesa@genome, names = aac_table$chr, start =  aac_table$nt1_pos - 1, end =  aac_table$nt1_pos + 1)
	second_triseq = BSgenome::getSeq(cesa@genome, names = aac_table$chr, start =  aac_table$nt2_pos - 1, end =  aac_table$nt2_pos + 1)
	third_triseq = BSgenome::getSeq(cesa@genome, names = aac_table$chr, start =  aac_table$nt3_pos - 1, end =  aac_table$nt3_pos + 1)

	# put together the three reference nucleotides of each codon (in coding order, not genomic order)
	coding_sequences = Biostrings::xscat(Biostrings::subseq(first_triseq, start = 2, width = 1),
	                                     Biostrings::subseq(second_triseq, start = 2, width = 1),
	                                     Biostrings::subseq(third_triseq, start = 2, width = 1))
	
	aac_table[, nt1_ref := as.character(Biostrings::subseq(coding_sequences, start = 1, width = 1))]
	aac_table[, nt2_ref := as.character(Biostrings::subseq(coding_sequences, start = 2, width = 1))]
	aac_table[, nt3_ref := as.character(Biostrings::subseq(coding_sequences, start = 3, width = 1))]
	
	
	# take complement of negative-strand coding sequences (but no need to reverse)
	coding_sequences[aac_table$strand == -1] = Biostrings::complement(coding_sequences[aac_table$strand == -1])
	aac_table[, coding_seq := as.character(coding_sequences)]
	
	# for each aa_mut, get a three-item list representing SNVs in first, second, third positions of the reference codon
	# that cause the same amino acid change as aa_alt (which will always include the ntchange that came from the MAF)
	possible_snvs = mapply(function(x, y) codon_snvs_to_aa[[x]][[y]], as.character(coding_sequences), aac_table$aa_alt, SIMPLIFY = F)
	
	
	## Sanity check
	# tmp = seqinr::aaa(as.character(Biostrings::translate(coding_sequences, no.init.codon = T)))
	# tmp[tmp == "Stp"] = "STOP"
	# identical(tmp, aac_table$aa_ref)
	# 
	
	nt1_snvs = lapply(possible_snvs, function(x) x[[1]])
	aac_table[, nt1_snvs := nt1_snvs]
	
	nt2_snvs = lapply(possible_snvs, function(x) x[[2]])
	aac_table[, nt2_snvs := nt2_snvs]
	
	nt3_snvs = lapply(possible_snvs, function(x) x[[3]])
	aac_table[, nt3_snvs := nt3_snvs]
	
	# create table of all snv_ids, then fix strand of alt alleles
	nt1_dt = aac_table[, .(chr, strand, pos = nt1_pos, ref = nt1_ref, alt = unlist(nt1_snvs)), by = "aac_id"]
	nt2_dt = aac_table[, .(chr, strand, pos = nt2_pos, ref = nt2_ref, alt = unlist(nt2_snvs)), by = "aac_id"]
	nt3_dt = aac_table[, .(chr, strand, pos = nt3_pos, ref = nt3_ref, alt = unlist(nt3_snvs)), by = "aac_id"]

	snv_table = rbind(nt1_dt, nt2_dt, nt3_dt)
	snv_table = snv_table[! is.na(alt)] # drop positions where there is no possible SNV that causes the amino acid change of interest
  snv_table[strand == -1, alt := seqinr::comp(alt, forceToLower = F)]
	snv_table[, strand := NULL]
	snv_table[, snv_id := paste0(chr, ':', pos, "_", ref, '>', alt)]
	
	# annotate amino acid mutations with all associated SNVs
	snvs_by_aa_mut = snv_table[, .(all_snv_ids = list(snv_id)), by = "aac_id"]
	aac_table[snvs_by_aa_mut, all_snv_ids := all_snv_ids, on = "aac_id"]
	
	# make SNV table records unique and annotate with all amino acid mutations
  snv_table = snv_table[, .(chr, pos, ref, alt, assoc_aa_mut = list(aac_id)), by = "snv_id"]

  snv_table = unique(snv_table, by = "snv_id")
  
  # add in noncoding MAF SNVs
  noncoding = MAF[Variant_Type == "SNV" & ! snv_id %in% unlist(aac_table$all_snv_ids)]
  noncoding = noncoding[, .(chr = Chromosome, pos = Start_Position, ref = Reference_Allele, alt = Tumor_Allele)]
  noncoding[, snv_id := paste0(chr, ':', pos, "_", ref, '>', alt)]
  setcolorder(noncoding, "snv_id")
  noncoding = unique(noncoding)
  noncoding[, assoc_aa_mut := list(as.list(NA_character_))] # to match coding records
  snv_table = rbind(snv_table, noncoding)
  
  
  
  # find the closest reference gene to each MAF record (ties are possible, including one record being found in multiple genes)
  snv_grs = GenomicRanges::makeGRangesFromDataFrame(snv_table, seqnames.field = "chr", start.field = "pos", end.field = "pos")
  nearest = as.data.table(GenomicRanges::distanceToNearest(snv_grs, gr_genes_in_data, select = "all"))
  
  # convert the "subjectHits" index returned by the distanceToNearest function to the corresponding gene name
  possible_genes = GenomicRanges::mcols(gr_genes_in_data)["names"][,1]
  nearest[, gene := nearest[, possible_genes[subjectHits]]]
  
  # remove all but one of multiple hits from one record to the same gene
  # this happens when a record overlaps multiple exons in the reference data from different transcripts
  nearest = nearest[! duplicated(nearest[, .(queryHits, gene)])] 
  
  # queryHits column gives SNV table row row, gene gives associated gene
  # some records have multiple matching genes; combine them into variable-length vector within each snv table row
  # also grab the distance from the first matching gene (distances are always equal on multiple hits, since we asked for the nearest gene)
  genes_by_snv_row = nearest[, .(genes = list(gene), dist = distance[1]), by = "queryHits"]
  snv_table[, genes := genes_by_snv_row$genes]
  
  # confirm that all dNdScv gene annotations appear in our own gene annotation
  setkey(aac_table, "aac_id")
  ## Should make this check but it needs to run faster
  # genes_from_dndscv = lapply(snv_table$assoc_aa_mut, function(x) aac_table[x, gene])
  # 
  # bad_anno = which(! sapply(1:length(genes_from_dndscv), function(x) all(genes_from_dndscv[[x]] %in% snv_table$genes[[x]])))
  # if(length(bad_anno) > 0) {
  #   warning(paste0("The following SNVs were not annotated consistently between dNdScv and cancereffectsizeR:\n", 
  #                  paste0(snv_table[bad_anno, snv_id], collapse = ", "), "\n",
  #                  "This could indicate a problem with custom RefCDS data. If you are using CES's built-in RefCDS,\n",
  #                  "please consider submitting a bug report."))
  # }
  # 
  # when distance >0, call the record intergenic
  snv_table[, intergenic := TRUE][genes_by_snv_row$dist == 0, intergenic := FALSE]

  #  get deconstructSigs-style trinuc context of each SNV ID
  genomic_context = BSgenome::getSeq(cesa@genome, snv_table$chr,start=snv_table$pos - 1,
                                     end=snv_table$pos + 1,
                                     as.character = TRUE)
  trinuc_mut_ids = paste0(genomic_context,":", snv_table$alt)
  snv_table[, trinuc_mut := trinuc_translator[trinuc_mut_ids, "deconstructSigs_format"]]
  
  
  # clean up aa table
  aac_table = aac_table[, .(aac_id, gene, strand, pid, aachange, aa_ref, aa_pos, aa_alt, next_to_splice, nt1_pos, nt2_pos, nt3_pos, coding_seq, all_snv_ids)]
  setcolorder(aac_table, c("aac_id", "gene", "aachange", "strand"))
	# If any trinucleotide mutation comes up NA--usually due to an ambiguous N in the genomic trinuc context--remove record from analysis
	bad_trinuc_context = which(is.na(snv_table$trinuc_mut))
	
	# in the future, variant type test removed (for now, only SNV annotations really matter)
	num_bad = length(bad_trinuc_context)
	if (num_bad > 0) {
	  bad_ids = snv_table[bad_trinuc_context, snv_id]
	  bad_trinuc_context_maf <- MAF[snv_id %in% bad_ids, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
	  MAF <- MAF[! snv_id %in% bad_ids,]
	  message(paste("Note:", num_bad, "MAF records were excluded from analysis because of ambiguous trinucleotide context (likely N's in the reference genome)."))
	  bad_trinuc_context_maf$Exclusion_Reason = "ambiguous_trinuc_context"
	  cesa@excluded = rbind(cesa@excluded, bad_trinuc_context_maf)
	  
	  # remove the bad record from SNV and aa tables
	  bad_aa = unlist(snv_table[bad_trinuc_context, assoc_aa_mut])
	  snv_table = snv_table[! bad_trinuc_context]
	  aac_table = aac_table[! aac_id %in% bad_aa]
	}
	
	
	
	# record which covered_regions granges cover each mutation
	snv_gr = GenomicRanges::makeGRangesFromDataFrame(snv_table, seqnames.field = "chr", start.field = "pos", end.field = "pos")
	
	# test each MAF locus against all coverage grs
	# this returns a data frame where rows match MAF rows, columns are T/F for each coverage gr
	is_covered = as.data.table(lapply(cesa@coverage, function(x) snv_gr %within% x))

	
	# get the names of coverage grs with coverage for each site (and add in genome, which covers every site)
	grs_with_coverage = apply(is_covered, 1, function(x) c(names(which(x == TRUE)), "genome"))
	
	# when all samples have same coverage apply "helpfully" returns a matrix, but we want a list
	if(! is(grs_with_coverage, "list")) {
	  grs_with_coverage = as.list(as.data.table(grs_with_coverage))
	}
	
	# Note that when exome+ coverage (see load_maf) is used, samples can have both "exome" and "exome+" associated with their mutations,
	# but the samples themselves are considered "exome+" (be careful not to double-count these if developing something new)
	snv_table[, covered_in := grs_with_coverage]
	
	# We're going to cheat and and say samples have coverage at aac sites if they have coverage on any of the three codon positions
	# in practice, there probably is coverage even if the coverage intervals disagree
	setkey(snv_table, "snv_id")
	
	# speed this up later, maybe
	aac_table[, covered_in := .(list(snv_table[all_snv_ids, unique(unlist(covered_in))])), by = "aac_id"]
	
	snv_maf = MAF[Variant_Type == "SNV"]
	indel_maf = MAF[Variant_Type != "SNV"]
	
	snv_maf = cbind(snv_maf, snv_table[snv_maf, .(genes, assoc_aa_mut, trinuc_mut, covered_in) , on = "snv_id"])
	MAF = rbind(snv_maf, indel_maf, fill = T)
	
	# set all non-SNV annotation fields to NA (pending future development)
	MAF[Variant_Type != "SNV", c(7:ncol(MAF)) := NA]
	MAF = MAF[order(Chromosome, Start_Position)]
	
	# drop annotmuts info since it's already been used here (and it takes up a lot of memory)
	lapply(cesa@dndscv_out_list, function(x) x$annotmuts = NULL)


	cesa@maf = MAF
	cesa@mutations = list(amino_acid_change = aac_table, snv = snv_table)
	return(cesa)
}
