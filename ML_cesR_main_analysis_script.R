# load in packages required for cancereffectsizeR ----
library(cancereffectsizeR)


# Load in mutation data ----

NCI_data <- read.delim(file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/TCGA.LUAD.mutect.0458c57f-316c-4a7c-9294-ccd11c97c2f9.DR-10.0.somatic.maf",stringsAsFactors = F,skip = 5,header = T)

YG_data <- read.delim(file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/mutationsTN_62_Lung_Adenocarcinoma_(Yale___MD_Anderson).maf",header = T,stringsAsFactors = F)

# Preprocess for analysis ----
# These functions will be automatically loaded into the package space,
# but let's source them into this working file.


# TODO: check all these preprocessing steps to make sure they are necessary and efficient

## convert NCI data to hg19 ----
source("R/hg_converter.R")

NCI_data <- hg_converter(chain = "~/Box Sync/ML_cancereffectsizeR_data/hg38ToHg19.over.chain",maf_to_convert =  NCI_data)

## merge the data frames along common headers ----
source("R/merging_NCI_and_local_MAF_files.R")

MAF_input <- merging_TCGA_and_local_MAFdata_function(NCI_data = NCI_data,Local_data = YG_data,check_for_same_tumor = T)

## add unique tumor names ----
# is this necessary anymore? Want to remove superfluous functions and code.
source("R/unique_tumor_addition.R")

MAF_input <- unique_tumor_addition_function(MAF.file = MAF_input,non.TCGA.characters.to.keep = "all")

## tumor allele adder ----
source("R/tumor_allele_adder.R")

MAF_input <- tumor_allele_adder(MAF = MAF_input)

## remove DNP, TNP ----
source("R/DNP_TNP_remover.R")

MAF_input <- DNP_TNP_remover(MAF = MAF_input,delete_recur = T)

MAF_input <- MAF_input[,c("Unique_patient_identifier","Chromosome","Start_Position","Reference_Allele","Tumor_allele")]
## save for easy recall -----
save(MAF_input, file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/LUAD_MAF_for_analysis.RData")


# End of preprocessing ----

# load("~/Box Sync/ML_cancereffectsizeR_data/LUAD/LUAD_MAF_for_analysis.RData")

# Start of analysis ----

# These values will be included in the main function header:
sample_ID_column = "Unique_patient_identifier"
chr_column = "Chromosome"
pos_column = "Start_Position"
ref_column = "Reference_Allele"
alt_column = "Tumor_allele"
covariate_file = "lung_pca"


# will use functions and sub-functions to complete each step
# have clean and consistent inputs/outputs at each step
# will eventually wrap this all in a single function, that takes preprocessed MAF and outputs selection results.

# 1. Find trinucleotide signature weights ----
## Use deconstructSigs (https://github.com/raerose01/deconstructSigs)
### Keep formatting consistent (i.e. A[A>T]G) with deconstructSigs output (not original cancereffectsizeR
###

source("R/deconstructSigs_input_preprocess.R")

MAF_input_deconstructSigs_preprocessed <- deconstructSigs_input_preprocess(MAF = MAF_input)

# find tumors with more than or = to 50 variants
substitution_counts <- table(MAF_input_deconstructSigs_preprocessed[,sample_ID_column])
tumors_with_50_or_more <- names(which(substitution_counts>=50))
tumors_with_less_than_50 <- setdiff(MAF_input_deconstructSigs_preprocessed[,sample_ID_column],tumors_with_50_or_more)

# This provides us with the number of each trinucleotide in each tumor
# Even the ones with less than 50 variants
# Useful for the nearest neighbor calculation
trinuc_breakdown_per_tumor <- deconstructSigs::mut.to.sigs.input(mut.ref =
                                                                   MAF_input_deconstructSigs_preprocessed,
                                                                 sample.id = sample_ID_column,
                                                                 chr = chr_column,
                                                                 pos = pos_column,
                                                                 ref = ref_column,
                                                                 alt = alt_column)

deconstructSigs_trinuc_string <- colnames(trinuc_breakdown_per_tumor)
save(deconstructSigs_trinuc_string,file = "data/deconstructSigs_trinuc_string.RData")

# df to store relative proportions of all trinucs in all tumors

trinuc_proportion_matrix <- matrix(data = NA,
                                   nrow = nrow(trinuc_breakdown_per_tumor),
                                   ncol = ncol(trinuc_breakdown_per_tumor))
rownames(trinuc_proportion_matrix) <- rownames(trinuc_breakdown_per_tumor)
colnames(trinuc_proportion_matrix) <- colnames(trinuc_breakdown_per_tumor)



for(tumor_name in 1:length(tumors_with_50_or_more)){
  signatures_output <- deconstructSigs::whichSignatures(tumor.ref = trinuc_breakdown_per_tumor,
                                                        signatures.ref = signatures.cosmic,
                                                        sample.id = tumors_with_50_or_more[tumor_name],
                                                        contexts.needed = TRUE,
                                                        tri.counts.method = 'exome2genome')

  trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] <- signatures_output$product/sum( signatures_output$product) #need it to sum to 1.

  # Not all trinuc weights in the cosmic dataset are nonzero for certain signatures
  # This leads to the rare occasion where a certain combination of signatues leads to ZERO
  # rate for particular trinucleotide contexts.
  # True rate is nonzero, as we do see those variants in those tumors, so renormalizing
  # the rates by adding the lowest nonzero rate to all the rates and renormalizing


  if(0 %in% trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],]){
    # finding the lowest nonzero rate
    lowest_rate <- min(trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],
                                 -which(trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],]==0)])

    # adding it to the rates
    trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] <- trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] + lowest_rate

    # renormalizing to 1
    trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] <-
      trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],] / sum(trinuc_proportion_matrix[tumors_with_50_or_more[tumor_name],])
  }

}

# 2. Find nearest neighbor to tumors with < 50 mutations, assign identical weights as neighbor ----

distance_matrix <- as.matrix(dist(trinuc_breakdown_per_tumor))

for(tumor_name in 1:length(tumors_with_less_than_50)){
  #find closest tumor that have over 50 mutations
  closest_tumor <- names(sort(distance_matrix[tumors_with_less_than_50[tumor_name],tumors_with_50_or_more]))[1]

  trinuc_proportion_matrix[tumors_with_less_than_50[tumor_name],] <- trinuc_proportion_matrix[closest_tumor,]

}

# now, trinuc_proportion_matrix has the proportion of all trinucs in every tumor.

# 3. Calculate mutation rates at gene-level (dndscv) ----

MAF_input <- cancereffectsizeR::wrongRef_dndscv_checker(mutations = MAF_input[,c(sample_ID_column,chr_column,pos_column,ref_column,alt_column)])

if(is.null(covariate_file)){
  data("covariates_hg19",package = "dndscv")
  genes_in_pca <- rownames(covs)
}else{
  this_cov_pca <- get(data(list=covariate_file, package="cancereffectsizeR"))
  genes_in_pca <- rownames(this_cov_pca$rotation)
}

path_to_library <- dir(.libPaths(),full.names=T)[grep(dir(.libPaths(),full.names=T),pattern="cancereffectsizeR")][1] # find the path to this package


data("RefCDS_TP53splice",package = "cancereffectsizeR")

dndscvout <- dndscv::dndscv(
  mutations = MAF_input,
  gene_list = genes_in_pca,
  cv = if(is.null(covariate_file)){ "hg19"}else{ this_cov_pca$rotation},
  refdb = paste(path_to_library,"/data/RefCDS_TP53splice.RData",sep=""))



data("refcds_hg19", package="dndscv")

RefCDS.our.genes <- RefCDS[which(sapply(RefCDS, function(x) x$gene_name) %in% dndscvout$genemuts$gene_name)]


if(dndscvout$nbreg$theta>1){
  mutrates <- ((dndscvout$genemuts$n_syn + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv))/sapply(RefCDS.our.genes, function(x) colSums(x$L)[1]))/length(unique(MAF_input[,sample_ID_column]))
}else{
  mutrates <- rep(NA,length(dndscvout$genemuts$exp_syn_cv))

  syn.sites <- sapply(RefCDS.our.genes, function(x) colSums(x$L)[1])

  for(i in 1:length(mutrates)){
    if( dndscvout$genemuts$exp_syn_cv[i] >  ((dndscvout$genemuts$n_syn[i] + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv[i])))){
      mutrates[i] <- (dndscvout$genemuts$exp_syn_cv[i]/syn.sites[i])/length(unique(MAF_input[,sample_ID_column]))
    }else{
      mutrates[i] <- (((dndscvout$genemuts$n_syn[i] + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv[i])))/syn.sites[i])/length(unique(MAF[,sample_ID_column]))
    }

  }
}
names(mutrates) <- dndscvout$genemuts$gene_name



# 4. Assign genes to MAF_input ----

MAF_ranges <- GenomicRanges::GRanges(seqnames = MAF_input[,chr_column], ranges = IRanges::IRanges(start=MAF_input[,pos_column],end = MAF_input[,pos_column]))

data("refcds_hg19", package="dndscv") # load in gr_genes data
gene_name_matches <- GenomicRanges::nearest(x = MAF_ranges, subject = gr_genes,select=c("all"))

MAF_input$Gene_name <- NA

# take care of single hits
single.choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))==1)])

MAF_input$Gene_name[single.choice] <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)[which(S4Vectors::queryHits(gene_name_matches) %in% single.choice)]]


multi.choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))>1)])

all.possible.names <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)]
query.spots <- S4Vectors::queryHits(gene_name_matches)

for(i in 1:length(multi.choice)){
  genes.for.this.choice <- all.possible.names[which(query.spots==multi.choice[i])]
  if(length(which( genes.for.this.choice %in% names(mutrates) ))>0){
    MAF_input$Gene_name[multi.choice[i]] <- genes.for.this.choice[which( genes.for.this.choice %in% names(mutrates) )[1]]
  }else{
    MAF_input$Gene_name[multi.choice[i]] <- "Indeterminate"
  }
}




# 5. For each substitution, calculate the gene- and tumor- and mutation-specific mutation rate----


load("data/gene_trinuc_comp.RData")
load("data/AA_mutation_list.RData")

# save.image(file = "~/Box Sync/ML_cancereffectsizeR_data/LUAD/before_step_5.RData")
# load("~/Box Sync/ML_cancereffectsizeR_data/LUAD/before_step_5.RData")

MAF_input <- MAF_input[which(MAF_input$Gene_name %in% names(mutrates)),]
MAF_input <- MAF_input[which(MAF_input$Reference_Allele %in% c("A","T","C","G") & MAF_input$Tumor_allele %in% c("A","T","C","G")),]

dndscvout_annotref <- dndscvout$annotmuts[which(dndscvout$annotmuts$ref %in% c("A","T","G","C") & dndscvout$annotmuts$mut %in% c("A","T","C","G")),]

# assign strand data
strand_data <- sapply(RefCDS, function(x) x$strand)
names(strand_data) <- sapply(RefCDS, function(x) x$gene_name)
strand_data[which(strand_data==-1)] <- "-"
strand_data[which(strand_data==1)] <- "+"

MAF_input$strand <- strand_data[MAF_input$Gene_name]

# strReverse <- function(x)
#   sapply(lapply(strsplit(x, NULL), rev), paste, collapse="") #found in ?strsplit

MAF_input$triseq <- as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste("chr",MAF_input$Chromosome,sep=""), strand=MAF_input$strand, start=MAF_input$Start_Position-1, end=MAF_input$Start_Position+1))

# MAF_input[which(MAF_input$strand=="-"),"triseq"] <- strReverse(MAF_input[which(MAF_input$strand=="-"),"triseq"])


MAF_input$unique_variant_ID <- paste(MAF_input$Chromosome,MAF_input$Start_Position, MAF_input$Tumor_allele)
dndscvout_annotref$unique_variant_ID <- paste(dndscvout_annotref$chr, dndscvout_annotref$pos, dndscvout_annotref$mut)

MAF_input$is_coding <- MAF_input$unique_variant_ID %in% dndscvout_annotref$unique_variant_ID[which(dndscvout_annotref$impact != "Essential_Splice")]

# Assign trinucleotide context data (for use with non-coding variants)
MAF_input$Tumor_allele_correct_strand <- MAF_input$Tumor_allele
MAF_input$Tumor_allele_correct_strand[which(MAF_input$strand=="-")] <- toupper(seqinr::comp(MAF_input$Tumor_allele_correct_strand[which(MAF_input$strand=="-")]))

load("data/trinuc_translator.RData")

MAF_input$trinuc_dcS <- trinuc_translator[paste(MAF_input$triseq,":",MAF_input$Tumor_allele_correct_strand,sep=""),"deconstructSigs_format"]

# Assign coding variant data

dndscv_coding_unique <- dndscvout_annotref[!duplicated(dndscvout_annotref$unique_variant_ID),]

rownames(dndscv_coding_unique) <- dndscv_coding_unique$unique_variant_ID
MAF_input$nuc_variant <- NA
MAF_input$coding_variant <- NA
MAF_input[,c("nuc_variant","coding_variant")] <- dndscv_coding_unique[MAF_input$unique_variant_ID,c("ntchange","aachange")]
MAF_input[which(MAF_input$coding_variant=="-"),"coding_variant"] <- NA
MAF_input[which(MAF_input$nuc_variant=="-"),"nuc_variant"] <- NA

MAF_input$nuc_position  <- as.numeric(gsub("\\D", "", MAF_input$nuc_variant))
MAF_input$codon_pos <- (MAF_input$nuc_position %% 3)
MAF_input$codon_pos[which(MAF_input$codon_pos==0)] <- 3

# TODO: problem with the way I generate the sequence to find all
# possible trinucleotides responsible for AA change observed. If codon lies within
# or adjacent to exon start/stop, then trinucleotides responsible for that
# variant is different than seen in neighbors considered




MAF_input$amino_acid_context <- as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste("chr",MAF_input$Chromosome,sep=""), strand=MAF_input$strand, start=MAF_input$Start_Position-3, end=MAF_input$Start_Position+3))

ref_cds_genes <- sapply(RefCDS, function(x) x$gene_name)
names(RefCDS) <- ref_cds_genes

# need to take CDS into account.

# issue: if the variant is at the end of the gene then a shorter sequence is produced.

# for(i in 1:nrow(MAF_input)){
#  if(MAF_input$is_coding[i]){
#   MAF_input$amino_acid_context[i] <- substr(as.character(RefCDS[[MAF_input$Gene_name[i]]]$seq_cds),MAF_input$nuc_position[i]-3,MAF_input$nuc_position[i]+3)
#   if(MAF_input$nuc_position[i] == RefCDS[[MAF_input$Gene_name[i]]]$CDS_length){
#     MAF_input$amino_acid_context[i] <- paste(MAF_input$amino_acid_context[i],substr(as.character(RefCDS[[MAF_input$Gene_name[i]]]$seq_cds1down),RefCDS[[MAF_input$Gene_name[i]]]$CDS_length,RefCDS[[MAF_input$Gene_name[i]]]$CDS_length),collapse = "",sep="")
#   }
#
#   if(MAF_input$nuc_position[i] <= 3){
#     MAF_input$amino_acid_context[i] <- paste(substr(as.character(RefCDS[[MAF_input$Gene_name[i]]]$seq_cds1up),1,1),substr(MAF_input$amino_acid_context[i],1,4) ,collapse = "",sep="")
#   }
#  }
#
# }




# MAF_input[which(MAF_input$strand=="-"),"amino_acid_context"] <- strReverse(MAF_input[which(MAF_input$strand=="-"),"amino_acid_context"])

MAF_input$amino_acid_context[which(MAF_input$codon_pos==1)] <- substr(MAF_input$amino_acid_context[which(MAF_input$codon_pos==1)],3,7)

MAF_input$amino_acid_context[which(MAF_input$codon_pos==2)] <- substr(MAF_input$amino_acid_context[which(MAF_input$codon_pos==2)],2,6)

MAF_input$amino_acid_context[which(MAF_input$codon_pos==3)] <- substr(MAF_input$amino_acid_context[which(MAF_input$codon_pos==3)],1,5)

MAF_input$next_to_splice <- F
for(i in 1:nrow(MAF_input)){
  if(any(abs(MAF_input$Start_Position[i] - RefCDS[[MAF_input$Gene_name[i]]]$intervals_cds) <= 3)){

    MAF_input$next_to_splice[i] <- T

  }
}


names(gene_trinuc_comp) <- sapply(RefCDS, function(x) x$gene_name)

MAF_input$unique_variant_ID_AA <- NA

MAF_input[which(MAF_input$is_coding==T),"unique_variant_ID_AA"] <- MAF_input[which(MAF_input$is_coding==T),"coding_variant"]
MAF_input[which(MAF_input$is_coding==F),"unique_variant_ID_AA"] <- MAF_input[which(MAF_input$is_coding==F),"unique_variant_ID"]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

MAF_input$coding_variant_AA_mut <-  substrRight(x = MAF_input$coding_variant,n=1)

data("AA_translations") # from cancereffectsizeR

AA_translations_unique <- AA_translations[!duplicated(AA_translations$AA_letter),]
rownames(AA_translations_unique) <- as.character(AA_translations_unique$AA_letter)

MAF_input$coding_variant_AA_mut <- as.character(AA_translations_unique[MAF_input$coding_variant_AA_mut,"AA_short"])

# Now time for the function that calculates gene- and tumor- and trinuc- specific rate.

source("R/mutation_rate_calc.R")

KRAS_muts <- mutation_rate_calc(MAF = MAF_input,gene = "KRAS",gene_mut_rate = mutrates,tumor_trinucs = trinuc_proportion_matrix)


# 6. Maximum-likelihood framework

tumors <- unique(MAF_input$Unique_patient_identifier)

source("R/selection_intensity_calc_function_ML.R")

optimize_gamma(MAF_input=MAF_input, all_tumors=tumors, gene="KRAS", variant="G12C", specific_mut_rates=KRAS_muts)

save.image("~/Box Sync/ML_cancereffectsizeR_data/LUAD/after_step_6.RData")


# test <- cancereffectsizeR::trinuc_profile_function_with_weights(input.MAF = MAF_input)






