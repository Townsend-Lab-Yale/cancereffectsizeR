#' Annotates the SNV MAF with gene information, keeping assignments consistent with dndscv when possible
#'
#'
#'
#' @param cesa CESAnalysis object
#' @export


annotate_gene_maf <- function(cesa) {
	MAF = cesa@maf
	bases = c("A", "C", "T", "G")

	# subset to SNVs
	MAF = MAF[MAF$Reference_Allele %in% bases & MAF$Tumor_Allele %in% bases,]
	RefCDS_our_genes = cesa@refcds_data
	dndscv_gene_names = names(cesa@mutrates_list[[1]])
	dndscv_out_list = cesa@dndscv_out_list
	MAF$Gene_name <- NA
	MAF$unsure_gene_name <- F

	MAF_ranges <- GenomicRanges::GRanges(seqnames = MAF$Chromosome, ranges = IRanges::IRanges(start=MAF$Start_Position,end = MAF$Start_Position))


	# first, find overlaps

	gene_name_overlaps <- GenomicRanges::findOverlaps(query = MAF_ranges,subject = gr_genes, type="any", select="all")

	# duplicate any substitutions that matched two genes
	overlaps <- as.matrix(gene_name_overlaps)

	matched_ol <- which(1:nrow(MAF) %in% overlaps[,1])
	unmatched_ol <- setdiff(1:nrow(MAF), matched_ol)

	MAF_matched <- MAF[matched_ol,]

	MAF_ranges <- GenomicRanges::GRanges(seqnames = MAF_matched$Chromosome, ranges = IRanges::IRanges(start=MAF_matched$Start_Position,end = MAF_matched$Start_Position))

	gene_name_overlaps <- as.matrix(GenomicRanges::findOverlaps(query = MAF_ranges,subject = gr_genes, type="any", select="all"))

	MAF_matched <- MAF_matched[gene_name_overlaps[,1],] #expand the multi-matches
	MAF_matched$Gene_name <- gr_genes$names[gene_name_overlaps[,2]]# assign the multi-matches


	MAF_unmatched <- MAF[unmatched_ol,]



	# MAF_expanded <- MAF[overlaps[,1],] #this gets rid of unmatched subs, though...

	# search out the unmatched remainder

	MAF_ranges <- GenomicRanges::GRanges(seqnames = MAF_unmatched$Chromosome, ranges = IRanges::IRanges(start=MAF_unmatched$Start_Position,end = MAF_unmatched$Start_Position))

	gene_name_matches <- GenomicRanges::nearest(x = MAF_ranges, subject = gr_genes,select=c("all"))


	# take care of single hits
	single_choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))==1)])

	MAF_unmatched$Gene_name[single_choice] <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)[which(S4Vectors::queryHits(gene_name_matches) %in% single_choice)]]


	# Then, assign rest to the closest gene
	multi_choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))>1)])
	all_possible_names <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)]
	query_spots <- S4Vectors::queryHits(gene_name_matches)




	for(i in 1:length(multi_choice)){
		# first, assign if the nearest happened to be within the GenomicRanges::findOverlaps() and applicable to one gene
		if(length(which(queryHits(gene_name_matches) == multi_choice[i])) == 1){

			MAF_unmatched[multi_choice[i],"Gene_name"] <- gr_genes$names[subjectHits(gene_name_matches)[which(queryHits(gene_name_matches) == multi_choice[i])]]

		}else{

			genes_for_this_choice <- all_possible_names[which(query_spots==multi_choice[i])]

			if(any( genes_for_this_choice %in% dndscv_gene_names, na.rm=TRUE )){
				MAF_unmatched$Gene_name[multi_choice[i]] <- genes_for_this_choice[which( genes_for_this_choice %in% dndscv_gene_names )[1]]
			}else{
				MAF_unmatched$Gene_name[multi_choice[i]] <- "Indeterminate"
			}
		}
	}

	MAF_unmatched$unsure_gene_name <- T
	MAF <- rbind(MAF_matched,MAF_unmatched)


	data("AA_mutation_list", package = "cancereffectsizeR")

	MAF <- MAF[which(MAF$Gene_name %in% dndscv_gene_names),]

	dndscvout_annotref <- NULL
	for(this_subset in 1:length(dndscv_out_list)){
	  dndscvout_annotref <- rbind(dndscvout_annotref,dndscv_out_list[[this_subset]]$annotmuts)
	}

	dndscvout_annotref <- dndscvout_annotref[which(dndscvout_annotref$ref %in% c("A","T","G","C") & dndscvout_annotref$mut %in% c("A","T","C","G")),]

	# assign strand data
	strand_data <- sapply(RefCDS_our_genes, function(x) x$strand)
	names(strand_data) <- sapply(RefCDS_our_genes, function(x) x$gene_name)
	strand_data[which(strand_data==-1)] <- "-"
	strand_data[which(strand_data==1)] <- "+"

	MAF$strand <- strand_data[MAF$Gene_name]

	MAF$triseq <- as.character(
	  BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
	                   paste("chr",MAF$Chromosome,sep=""),
	                   strand=MAF$strand,
	                   start=MAF$Start_Position-1,
	                   end=MAF$Start_Position+1))

	# If trinucleotide context yields N (or other non-ACTG character), remove record from analysis
	bad_trinuc_context <- grepl('[^ACTG]', MAF$triseq)
	if (any(bad_trinuc_context)) {
		bad_trinuc_context_maf <- MAF[bad_trinuc_context,]
		MAF <- MAF[!bad_trinuc_context,]
		num_bad_records = nrow(bad_trinuc_context_maf) 
		message(paste("Note:", num_bad_records, "MAF records were excluded from analysis because the reference genome has N's (non-specific bases) in their trinucleotide context."))
		
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
	MAF$Tumor_allele_correct_strand <- MAF$Tumor_Allele
	MAF$Tumor_allele_correct_strand[which(MAF$strand=="-")] <-
	  toupper(seqinr::comp(MAF$Tumor_allele_correct_strand[which(MAF$strand=="-")]))

	data("trinuc_translator", package="cancereffectsizeR")

	MAF$trinuc_dcS <- trinuc_translator[paste(MAF$triseq,
	                                          ":",
	                                          MAF$Tumor_allele_correct_strand,
	                                          sep=""),"deconstructSigs_format"]

	# Assign coding variant data

	dndscv_coding_unique <- dndscvout_annotref[!duplicated(dndscvout_annotref$unique_variant_ID),]

	rownames(dndscv_coding_unique) <- dndscv_coding_unique$unique_variant_ID
	MAF$nuc_variant <- NA
	MAF$coding_variant <- NA
	MAF[,c("nuc_variant","coding_variant")] <- dndscv_coding_unique[MAF$unique_variant_ID,c("ntchange","aachange")]
	MAF[which(MAF$coding_variant=="-"),"coding_variant"] <- NA
	MAF[which(MAF$nuc_variant=="-"),"nuc_variant"] <- NA

	MAF$nuc_position  <- as.numeric(gsub("\\D", "", MAF$nuc_variant))
	MAF$codon_pos <- (MAF$nuc_position %% 3)
	MAF$codon_pos[which(MAF$codon_pos==0)] <- 3


	MAF$amino_acid_context <- as.character(
	  BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste(
	    "chr",MAF$Chromosome,sep=""),
	    strand=MAF$strand, start=MAF$Start_Position-3, end=MAF$Start_Position+3))

	ref_cds_genes <- sapply(RefCDS_our_genes, function(x) x$gene_name)
	names(RefCDS_our_genes) <- ref_cds_genes

	MAF$amino_acid_context[which(MAF$codon_pos==1)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==1)],3,7)

	MAF$amino_acid_context[which(MAF$codon_pos==2)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==2)],2,6)

	MAF$amino_acid_context[which(MAF$codon_pos==3)] <-
	  substr(MAF$amino_acid_context[which(MAF$codon_pos==3)],1,5)

	MAF$next_to_splice <- F
	for(i in 1:nrow(MAF)){
	  if(any(abs(MAF$Start_Position[i] - RefCDS_our_genes[[MAF$Gene_name[i]]]$intervals_cds) <= 3)){

	    MAF$next_to_splice[i] <- T

	  }
	}



	MAF$unique_variant_ID_AA <- NA

	MAF[which(MAF$is_coding==T),"unique_variant_ID_AA"] <- MAF[which(MAF$is_coding==T),"coding_variant"]
	MAF[which(MAF$is_coding==F),"unique_variant_ID_AA"] <- MAF[which(MAF$is_coding==F),"unique_variant_ID"]

	substrRight <- function(x, n){
	  substr(x, nchar(x)-n+1, nchar(x))
	}

	MAF$coding_variant_AA_mut <-  substrRight(x = MAF$coding_variant,n=1)

	data("AA_translations", package="cancereffectsizeR") # from cancereffectsizeR

	AA_translations_unique <- AA_translations[!duplicated(AA_translations$AA_letter),]
	rownames(AA_translations_unique) <- as.character(AA_translations_unique$AA_letter)

	MAF$coding_variant_AA_mut <- as.character(AA_translations_unique[MAF$coding_variant_AA_mut,"AA_short"])

	cesa@annotated.snv.maf = MAFdf(MAF)
	cesa@refcds_data = RefCDS_our_genes # has been updated by this script
	return(cesa)
}
