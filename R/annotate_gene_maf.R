#' Annotates the SNV MAF with gene information, keeping assignments consistent with dndscv when possible
#'
#'
#'
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
	

	MAF$triseq <- as.character(
	  BSgenome::getSeq(cesa@genome,
	                   MAF$Chromosome,
	                   strand=MAF$strand,
	                   start=MAF$Start_Position-1,
	                   end=MAF$Start_Position+1))

	# If trinucleotide context yields N (or other non-ACTG character), remove record from analysis
	bad_trinuc_context <- grepl('[^ACTG]', MAF$triseq)
	if (any(bad_trinuc_context)) {
		bad_trinuc_context_maf <- MAF[bad_trinuc_context, 1:5]
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
	
	MAF$Tumor_allele_correct_strand[which(MAF$strand=="-")] <-
	  toupper(seqinr::comp(MAF$Tumor_allele_correct_strand[which(MAF$strand=="-")]))


	MAF$trinuc_dcS <- trinuc_translator[paste(MAF$triseq,
	                                          ":",
	                                          MAF$Tumor_allele_correct_strand,
	                                          sep=""),"deconstructSigs_format"]

	# Assign coding variant data
	dndscv_coding_unique <- dndscvout_annotref[!duplicated(unique_variant_ID)]

	rownames(dndscv_coding_unique) <- dndscv_coding_unique$unique_variant_ID
	MAF$nuc_variant <- character(nrow(MAF))
	MAF$coding_variant <- character(nrow(MAF))
	MAF[,c("nuc_variant","coding_variant")] <- dndscv_coding_unique[MAF, .(as.character(ntchange),as.character(aachange)), on = "unique_variant_ID"]
	MAF[which(MAF$coding_variant=="-"),"coding_variant"] <- NA
	MAF[which(MAF$nuc_variant=="-"),"nuc_variant"] <- NA

	MAF$nuc_position  <- as.numeric(gsub("\\D", "", MAF$nuc_variant))
	MAF$codon_pos <- (MAF$nuc_position %% 3)
	MAF$codon_pos[which(MAF$codon_pos==0)] <- 3

	MAF$amino_acid_context <- as.character(
	  BSgenome::getSeq(cesa@genome, MAF$Chromosome, strand=MAF$strand, 
	                   start=MAF$Start_Position-3, end=MAF$Start_Position+3))

	MAF$amino_acid_context[which(MAF$codon_pos==1)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==1)],3,7)

	MAF$amino_acid_context[which(MAF$codon_pos==2)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==2)],2,6)

	MAF$amino_acid_context[which(MAF$codon_pos==3)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==3)],1,5)

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


	# drop annotmuts info since it's already been used here (and it takes up a lot of memory)
	lapply(cesa@dndscv_out_list, function(x) x$annotmuts = NULL)

	cesa@annotated.snv.maf = MAF
	return(cesa)
}
