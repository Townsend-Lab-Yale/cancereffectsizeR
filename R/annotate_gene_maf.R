#' Annotates the SNV MAF with gene information, keeping assignments consistent with dndscv when possible
#'
#'
#' @import 
#' @param cesa CESAnalysis object
#' @export
annotate_gene_maf <- function(cesa) {
  RefCDS = get_genome_data(cesa, "RefCDS")
  gr_genes = get_genome_data(cesa, "gr_genes")
  
  MAF = cesa@maf
	bases = c("A", "C", "T", "G")
  
	# subset to SNVs
	MAF = MAF[MAF$Reference_Allele %in% bases & MAF$Tumor_Allele %in% bases,]
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
	dndscvout_annotref <- dndscvout_annotref[ref %in% bases & mut %in% bases]

	# assign strand data
	strand_data = rbindlist(lapply(RefCDS, function(x) return(list(Gene_name = x$gene_name, strand = x$strand))))
	strand_data[, strand := as.character(strand)]
	strand_data[strand == "-1", strand := "-"]
	strand_data[strand == "1", strand := "+"]
	MAF[, strand := strand_data[.SD, strand, on = .(Gene_name)]]
	

	MAF$genomic_ref_template_triseq <- BSgenome::getSeq(cesa@genome,
	                   MAF$Chromosome,
	                   strand=MAF$strand,
	                   start=MAF$Start_Position-1,
	                   end=MAF$Start_Position+1, as.character = TRUE)

	# If trinucleotide context yields N (or other non-ACTG character), remove record from analysis
	bad_trinuc_context <- grepl('[^ACTG]', MAF$genomic_ref_template_triseq)
	if (any(bad_trinuc_context)) {
		bad_trinuc_context_maf <- MAF[bad_trinuc_context, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele)]
		MAF <- MAF[!bad_trinuc_context,]
		num_bad_records = nrow(bad_trinuc_context_maf) 
		message(paste("Note:", num_bad_records, "MAF records were excluded from analysis because the reference genome has N's (non-specific bases) in their trinucleotide context."))
		bad_trinuc_context_maf$Exclusion_Reason = "ambiguous_trinuc_context"
		cesa@excluded = rbind(cesa@excluded, bad_trinuc_context_maf)
	}

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
  
	# Assign trinucleotide context data (for use with non-coding variants)
	MAF[, Tumor_allele_correct_strand := Tumor_Allele]
	MAF[strand == "-", Tumor_allele_correct_strand := seqinr::comp(Tumor_allele_correct_strand, forceToLower = F)]
	MAF$trinuc_dcS <- trinuc_translator[paste0(MAF$genomic_ref_template_triseq,":",MAF$Tumor_allele_correct_strand), "deconstructSigs_format"]
	

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
	upstream = DNAStringSet(lapply(RefCDS[genes], function(x) x$seq_cds1up))
	inframe = DNAStringSet(lapply(RefCDS[genes], function(x) x$seq_cds))
	downstream = DNAStringSet(lapply(RefCDS[genes], function(x) x$seq_cds1down))
	
	# get the first, second, and third nucleotides of the given codon, and the genomic trinuc context of each
	# will shift from the given codon position to capture all three nucleotides in the codon
	start = coding_maf$nuc_position - coding_maf$codon_pos + 1
	upstream = Biostrings::subseq(upstream, start = start, width = 3)
	ref_codon = Biostrings::subseq(inframe, start = start, width = 3)
	downstream = Biostrings::subseq(downstream, start = start, width = 3)
	
	
	
	# handle non-splice
  not_splice = coding_maf[, next_to_splice == F]
  amino_acid_context = as.character(xscat(subseq(upstream[not_splice], start = 1, width = 1), 
                                          ref_codon[not_splice], 
                                          subseq(downstream[not_splice], start = 3, width = 1)))
  aa_mut_nonsplice = coding_maf[next_to_splice == F, coding_variant_AA_mut]
  equivalent_muts_nonsplice = lapply(1:length(amino_acid_context), function(x) unname(unlist(AA_mutation_list[[amino_acid_context[x]]][aa_mut_nonsplice[x]])))
  
  # handle splice
  splice = coding_maf[, next_to_splice == T]
  
  
  # for splice sites, calculate all possible SNV mutations (with trinuc context) that would match the actual amino acid change
  # this could also be pre-computed, but that hasn't happened yet
  # replace each position in the codon with all bases to get all codons that can results from an SNV mutation
  first_tri = xscat(subseq(upstream[splice], start = 1, width = 1), 
                    subseq(ref_codon[splice], start = 1, width = 1), 
                    subseq(downstream[splice], start = 1, width = 1))
  second_tri = xscat(subseq(upstream[splice], start = 2, width = 1),
                     subseq(ref_codon[splice], start = 2, width = 1),
                     subseq(downstream[splice], start = 2, width = 1))
  third_tri = xscat(subseq(upstream[splice], start = 3, width = 1),
                    subseq(ref_codon[splice], start = 3, width = 1),
                    subseq(downstream[splice], start = 3, width = 1))
	
  coding_splice_maf = coding_maf[next_to_splice == T, ]
  coding_splice_ref_codons = ref_codon[splice]
  
  equivalent_muts_splice = list()
  for(i in 1:nrow(coding_splice_maf)) {
    trinuc_contexts = DNAStringSet(list(first_tri[[i]], second_tri[[i]], third_tri[[i]]))
    codon = coding_splice_ref_codons[[i]]
    all_mutated_codons = list()
    j = 1
    # pre-calculate later
    bases = c("A", "T", "C", "G")
    for (nt in 1:4) {
      for(pos in 1:3) {
        all_mutated_codons[[j]] = replaceAt(codon, IRanges(pos, pos), bases[nt])
        j = j + 1
      }
    }
    # drop the non-mutated codon and duplicate entries
    all_mutated_codons = unique(sapply(all_mutated_codons, as.character))
    all_mutated_codons = all_mutated_codons[all_mutated_codons != as.character(codon)]
    
    # get translated amino acids for each SNV-mutated codon, and select those that match the actual amino acid mutation
    poss_amino_acids = sapply(all_mutated_codons, function(x) AA_translations[AA_translations$Nucs == x, "AA_short"])
    poss_mut = all_mutated_codons[which(poss_amino_acids == coding_splice_maf$coding_variant_AA_mut[i])]
    
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

	MAF[is_coding == FALSE, equivalent_aa_muts :=.(as.list(trinuc_dcS))]
	MAF[is_coding == TRUE & next_to_splice == FALSE, equivalent_aa_muts := equivalent_muts_nonsplice]
	MAF[is_coding == TRUE & next_to_splice == TRUE, equivalent_aa_muts := equivalent_muts_splice]

	# drop annotmuts info since it's already been used here (and it takes up a lot of memory)
	lapply(cesa@dndscv_out_list, function(x) x$annotmuts = NULL)


	## For comparing new/old annotated MAFs to ensure same selection results (see also mutation_rate_calc.R)
	MAF$amino_acid_context <- as.character(
	  BSgenome::getSeq(cesa@genome, MAF$Chromosome, strand=MAF$strand, 
	                   start=MAF$Start_Position-3, end=MAF$Start_Position+3))

	MAF$amino_acid_context[which(MAF$codon_pos==1)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==1)],3,7)

	MAF$amino_acid_context[which(MAF$codon_pos==2)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==2)],2,6)

	MAF$amino_acid_context[which(MAF$codon_pos==3)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==3)],1,5)


	cesa@annotated.snv.maf = MAF
	return(cesa)
}
