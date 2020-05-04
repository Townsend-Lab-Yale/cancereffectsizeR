#' annotate_variants
#' 
#' Annotates CESAnaysis MAF data with gene information, keeping assignments consistent with dndscv when possible
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
	dndscv_gene_names = names(cesa@mutrates_list[[1]])
	dndscv_out_list = cesa@dndscv_out_list

	# Select just the reference genes that are in the data output from dNdScv
	is_in_dndscv = (GenomicRanges::mcols(gr_genes)["names"][,1] %in% dndscv_gene_names)
	gr_genes_in_data = gr_genes[is_in_dndscv]
	
	# find the closest reference gene to each MAF record (ties are possible, including one record being found in multiple genes)
	gr_maf_data = GenomicRanges::GRanges(seqnames = MAF$Chromosome, ranges = IRanges::IRanges(start = MAF$Start_Position, end = MAF$Start_Position))
	nearest = as.data.table(GenomicRanges::distanceToNearest(gr_maf_data, gr_genes_in_data, select = "all"))
	
	# convert the "subjectHits" index returned by the distanceToNearest function to the corresponding gene name
	possible_genes = GenomicRanges::mcols(gr_genes_in_data)["names"][,1]
	nearest[,Gene_name := nearest[, possible_genes[subjectHits]]]
	
	# remove all but one of multiple hits from one record to the same gene
	# this happens, for example, when a record overlaps multiple exons in the reference data
	nearest = nearest[! duplicated(nearest[, .(queryHits, Gene_name)])] 
	
	# Note that some queryHits (corresponding to MAF row numbers) may have more than 1 hit, if 
	# a record fits in multiple genes (or, rarely, an exact tie for distance to nearest gene, if record isn't in any gene).
	# In these cases, a duplicate row of MAF is created for each additional gene.
	MAF = MAF[nearest$queryHits]
	MAF[, Gene_name := nearest$Gene_name]
	MAF[, unsure_gene_name := TRUE][nearest$distance == 0, unsure_gene_name := FALSE]

	# get mutation annotations from dNdScv and subset to SNVs
	dndscvout_annotref <- rbindlist(lapply(dndscv_out_list, function(x) x$annotmuts))
	dndscvout_annotref <- dndscvout_annotref[ref %in% c("A", "C", "G", "T") & mut %in% c("A", "C", "G", "T")]


	MAF$unique_variant_ID <- paste(
	  MAF$Gene_name,
	  MAF$Chromosome,
	  MAF$Start_Position,
	  MAF$Tumor_Allele)

	dndscvout_annotref$unique_variant_ID <-
	  paste(dndscvout_annotref$gene,
	        dndscvout_annotref$chr,
	        dndscvout_annotref$pos,
	        dndscvout_annotref$mut)

	MAF$is_coding <-
	  MAF$unique_variant_ID %in%
	  dndscvout_annotref$unique_variant_ID[which(dndscvout_annotref$impact != "Essential_Splice")]
  

	# Assign coding variant data
	dndscv_coding_unique <- dndscvout_annotref[!duplicated(unique_variant_ID)]

	rownames(dndscv_coding_unique) <- dndscv_coding_unique$unique_variant_ID
	MAF[,c("nuc_variant","coding_variant")] <- dndscv_coding_unique[MAF, .(as.character(ntchange),as.character(aachange)), on = "unique_variant_ID"]
	MAF[coding_variant == "-", coding_variant := NA]
	MAF[nuc_variant=="-", nuc_variant := NA]

	# get CDS intervals for every gene; potential splice sites are all the start/end positions of these
	splice_sites = rbindlist(lapply(RefCDS, function(x) return(list(Gene_name = x$gene_name, pos = x$intervals_cds))))
	comb = merge.data.table(MAF[,.(Gene_name, Start_Position)], splice_sites, allow.cartesian = T)
	comb[, diff:= abs(Start_Position - pos)]
	comb = comb[, any(diff <= 3), by = .(Gene_name, Start_Position)]
	MAF[, next_to_splice := comb[.SD, V1, on = .(Gene_name, Start_Position)]]
	
	
	MAF$unique_variant_ID_AA = character(nrow(MAF))
	
	MAF[is_coding==TRUE]$unique_variant_ID_AA = MAF[is_coding==TRUE]$coding_variant
	MAF[is_coding==FALSE]$unique_variant_ID_AA = MAF[is_coding==FALSE]$unique_variant_ID
	
	
	substrRight <- function(x, n){
	  substr(x, nchar(x)-n+1, nchar(x))
	}
	
	MAF$coding_variant_AA_mut <-  substrRight(x = MAF$coding_variant,n=1)
	AA_translations_unique <- AA_translations[!duplicated(AA_translations$AA_letter),]
	rownames(AA_translations_unique) <- as.character(AA_translations_unique$AA_letter)
	MAF$coding_variant_AA_mut <- as.character(AA_translations_unique[MAF$coding_variant_AA_mut,"AA_short"])
	
	
	MAF$nuc_position  <- as.numeric(gsub("\\D", "", MAF$nuc_variant))
	MAF[, codon_pos := nuc_position %% 3]
	MAF[codon_pos==0, codon_pos := 3] # "codon 0" obtained from  (nuc_position % 3) is actually codon 3

	# calculate genomic trinuc contexts for each base in the codon of each coding variant
	# transcript reference data is used so that variants near splice sites get handled appropriately
	coding_maf = MAF[is_coding == TRUE]
	genes = coding_maf[,Gene_name]
	upstream = Biostrings::DNAStringSet(lapply(RefCDS[genes], function(x) x$seq_cds1up))
	inframe = Biostrings::DNAStringSet(lapply(RefCDS[genes], function(x) x$seq_cds))
	downstream = Biostrings::DNAStringSet(lapply(RefCDS[genes], function(x) x$seq_cds1down))
	
	# get the first, second, and third nucleotides of the given codons, and the genomic trinuc context of each
	# will shift from the given codon position to capture all three nucleotides in the codon
	start = coding_maf$nuc_position - coding_maf$codon_pos + 1
	upstream = Biostrings::subseq(upstream, start = start, width = 3)
	ref_codon = Biostrings::subseq(inframe, start = start, width = 3)
	downstream = Biostrings::subseq(downstream, start = start, width = 3)
	
	# trinuc context for annotation (context for the mutated position, not the codon)
	transcript_trinuc = Biostrings::xscat(Biostrings::subseq(upstream, start = coding_maf$codon_pos, width = 1),
	                                      Biostrings::subseq(ref_codon, start = coding_maf$codon_pos, width = 1),
	                                      Biostrings::subseq(downstream, start = coding_maf$codon_pos, width = 1))
	strands = sapply(RefCDS[genes], function(x) x$strand)
	
	# take reverse complement in order to get 5' to 3' orientation for negative-strand genes
	genomic_trinuc = transcript_trinuc
	genomic_trinuc[strands == -1] = Biostrings::reverseComplement(transcript_trinuc[strands == -1])
	trinuc_dcS_coding <- trinuc_translator[paste0(genomic_trinuc,":",coding_maf$Tumor_Allele), "deconstructSigs_format"]
	
	# handle coding non-splice
  not_splice = coding_maf[, next_to_splice == F]
  amino_acid_context = as.character(Biostrings::xscat(Biostrings::subseq(upstream[not_splice], start = 1, width = 1), 
                                                                         ref_codon[not_splice], 
                                                      Biostrings::subseq(downstream[not_splice], start = 3, width = 1)))
  aa_mut_nonsplice = coding_maf[next_to_splice == F, coding_variant_AA_mut]
  equivalent_muts_nonsplice = lapply(1:length(amino_acid_context), function(x) unname(unlist(AA_mutation_list[[amino_acid_context[x]]][aa_mut_nonsplice[x]])))
  
  # handle coding splice
  splice = coding_maf[, next_to_splice == T]
  
  # for splice sites, calculate all possible SNV mutations (with trinuc context) that would match the actual amino acid change
  # replace each position in the codon with all bases to get all codons that can results from an SNV mutation
  first_tri = Biostrings::xscat(Biostrings::subseq(upstream[splice], start = 1, width = 1), 
                    Biostrings::subseq(ref_codon[splice], start = 1, width = 1), 
                    Biostrings::subseq(downstream[splice], start = 1, width = 1))
  second_tri = Biostrings::xscat(Biostrings::subseq(upstream[splice], start = 2, width = 1),
                     Biostrings::subseq(ref_codon[splice], start = 2, width = 1),
                     Biostrings::subseq(downstream[splice], start = 2, width = 1))
  third_tri = Biostrings::xscat(Biostrings::subseq(upstream[splice], start = 3, width = 1),
                    Biostrings::subseq(ref_codon[splice], start = 3, width = 1),
                    Biostrings::subseq(downstream[splice], start = 3, width = 1))
	
  coding_splice_maf = coding_maf[next_to_splice == T, ]
  coding_splice_ref_codons = ref_codon[splice]
  
  equivalent_muts_splice = list()
  for(i in 1:nrow(coding_splice_maf)) {
    trinuc_contexts = Biostrings::DNAStringSet(list(first_tri[[i]], second_tri[[i]], third_tri[[i]]))
    codon = coding_splice_ref_codons[[i]] # codon is a Biostring here

    ## dictionary gives all codons that can be made from one point mutation to the original codon
    ## keys = original codon; values = list where each name is a new codon, value is the resulting amino acid 
    poss_mut = codon_point_mutation_dict[[as.character(codon)]]
    mut_is_match = which(poss_mut == coding_splice_maf$coding_variant_AA_mut[i])
    poss_mut = names(poss_mut)[mut_is_match]
    
	  # for each possible mutated codon, saves its trinuc-context representation in deconstructSigs format
	  trinucs_for_rate = character()
	  codon = strsplit(as.character(codon), "")[[1]]
	  for (mut_codon in poss_mut) {
	    mut_codon = strsplit(mut_codon, "")[[1]]
	    index = which(mut_codon != codon)
	    ref_tri = as.character(trinuc_contexts[[index]]) # get the genomic (not transcript) reference trinucleotide context for mutated site
	    mutated_base = mut_codon[index]
	    ds_format = trinuc_translator[paste(ref_tri, mutated_base, sep=":"), "deconstructSigs_format"]
	    trinucs_for_rate = c(trinucs_for_rate, ds_format)
	  }
	  equivalent_muts_splice[[i]] = trinucs_for_rate
  }

  # for non-coding mutations, query the genome to get trinuc-context point mutations in deconstructSigs format (e.g., "G[T>C]T")
  noncoding_maf = MAF[is_coding == F]
  triseq <- BSgenome::getSeq(cesa@genome,
                             noncoding_maf$Chromosome,
                             start=noncoding_maf$Start_Position-1,
                             end=noncoding_maf$Start_Position+1,
                             as.character = TRUE)
  
  trinuc_dcS_noncoding <- trinuc_translator[paste0(triseq,":",noncoding_maf$Tumor_Allele), "deconstructSigs_format"]
  
  MAF[is_coding == TRUE, trinuc_dcS := trinuc_dcS_coding]
  MAF[is_coding == FALSE, trinuc_dcS := trinuc_dcS_noncoding]
  
  # add all to the annotated MAF
	MAF[is_coding == FALSE, equivalent_aa_muts :=.(as.list(trinuc_dcS_noncoding))]
	MAF[is_coding == TRUE & next_to_splice == FALSE, equivalent_aa_muts := equivalent_muts_nonsplice]
	MAF[is_coding == TRUE & next_to_splice == TRUE, equivalent_aa_muts := equivalent_muts_splice]

	# If any trinucleotide mutation comes up NA--usually due to an ambiguous N in the genomic trinuc context--remove record from analysis
	bad_trinuc_context = sapply(MAF$equivalent_aa_muts, function(x) any(is.na(x)))
	# in the future, variant type test removed (for now, only SNV annotations really matter)
	bad_trinuc_context = bad_trinuc_context & MAF$Variant_Type == "SNV"
	num_bad = sum(bad_trinuc_context)
	if (num_bad > 0) {
	  bad_trinuc_context_maf <- MAF[bad_trinuc_context, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
	  MAF <- MAF[!bad_trinuc_context,]
	  message(paste("Note:", num_bad, "MAF records were excluded from analysis because of ambiguous trinucleotide context (likely N's in the reference genome)."))
	  bad_trinuc_context_maf$Exclusion_Reason = "ambiguous_trinuc_context"
	  cesa@excluded = rbind(cesa@excluded, bad_trinuc_context_maf)
	}
	
	
	# record which covered_regions granges cover each mutation
	maf_gr = GenomicRanges::makeGRangesFromDataFrame(MAF, seqnames.field = "Chromosome", start.field = "Start_Position",  end.field = "Start_Position")
	
	# test each MAF locus against all coverage grs
	# this returns a data frame where rows match MAF rows, columns are T/F for each coverage gr
	is_covered = as.data.table(lapply(cesa@coverage, function(x) maf_gr %within% x))

	
	# get the names of coverage grs with coverage for each site (and add in genome, which covers every site)
	grs_with_coverage = apply(is_covered, 1, function(x) c(names(which(x == TRUE)), "genome"))
	
	# when all samples have same coverage apply "helpfully" returns a matrix, but we want a list
	if(! is(grs_with_coverage, "list")) {
	  grs_with_coverage = as.list(as.data.table(grs_with_coverage))
	}
	
	# Note that when exome+ coverage (see load_maf) is used, samples can have both "exome" and "exome+" associated with their mutations,
	# but the samples thmeslves are considered "exome+" (be careful not to double-count these if developing something new)
	MAF[,covered_in := grs_with_coverage]
	
	# set all non-SNV annotation fields to NA (pending future development)
	MAF[Variant_Type != "SNV", c(7:ncol(MAF)) := NA]
	
	# drop annotmuts info since it's already been used here (and it takes up a lot of memory)
	lapply(cesa@dndscv_out_list, function(x) x$annotmuts = NULL)


	cesa@maf = MAF
	return(cesa)
}
