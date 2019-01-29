#' SNV effect size calculation function
#'
#' Function that determines mutation rate at the gene- and trinucleotide- level
#' and then calculates cancer effect size.
#' Works best on datasets with a large number of tumors, each with a large
#' number of SNV (The \code{deconstructSigs} package, used here to calculate
#' trinucleotide signatues, suggests at least 50 SNV per tumor used in the
#' signature calculation). Requres data to be in hg coordinates, you can use
#' the \code{hg_converter()} function provided with this package to convert
#' from one coordinate system to another.
#'
#' @import deconstructSigs
#' @import dndscv
#'
#' @param MAF_file MAF file with substitution data
#' @param covariate_file Either NULL, which runs dndscv without covariates, "dndscv_default" which uses the \code{dndscv}
#' default covariates, or one of these:  "bladder_pca"  "breast_pca"
#' "cesc_pca" "colon_pca" "esca_pca" "gbm_pca" "hnsc_pca" "kidney_pca" "lihc_pca" "lung_pca" "ov_pca" "pancreas_pca" "prostate_pca" "rectum_pca" "skin_pca"  "stomach_pca"  "thca_pca" "ucec_pca"
#' @param output_from_mainMAF the output from a previous selection run with all the
#' mutation and selection data, useful for bootstrapping because it speeds up the
#' analysis (no longer need to look up mutational context, etc.)
#' @param genes_for_effect_size_analysis genes to calculate effect sizes within.
#' @param sample_ID_column column in MAF with sample ID data
#' @param ref_column column in MAF with reference allele data
#' @param alt_column column in MAF with alternative allele data
#' @param pos_column column in MAF with chromosome nucleotide location data
#' @param chr_column column in MAF with chromosome data
#' @param custom_genome T if you are using a genome other than hg19
#' @param custom_location_RefCDS location of RefCDS location of custom genome (see dndscv user guide for more information http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/buildref.html
#' @param custom_BSgenome BSgenome you have installed and loaded prior to runnig the function
#'
#'
#' @export




#### all in R

# rm(list=ls())

effect_size_SNV <- function(MAF_file,
                            covariate_file=NULL,
                            output_from_mainMAF=NULL,
                            sample_ID_column="Unique_patient_identifier",
                            chr_column = "Chromosome",
                            pos_column = "Start_Position",
                            ref_column = "Reference_Allele",
                            alt_column = "Tumor_allele",
                            genes_for_effect_size_analysis="all",
                            custom_genome=F,custom_location_RefCDS=NULL,custom_BSgenome=NULL){

  # library("seqinr")
  # library("Biostrings")
  # library("MASS")
  # library("GenomicRanges")
  # library("dndscv")
  # library("reshape2")
  # library("BSgenome.Hsapiens.UCSC.hg19")
  # source("SNV_effect_size_functions")
  # dndscv column names sampleID, chr, pos, ref, alt


  # source("../R/dndscv_wrongRef_checker.R")
  # source("../R/all_in_R_trinucleotide_profile.R")
  # source("../R/selection_intensity_calc_dndscv.R")

  # load("../R/all_gene_trinuc_data.RData")

  # load in MAF

  # MAF <- get(load(MAF_file))

  if(custom_genome){
    if(is.null(custom_location_RefCDS) | is.null(custom_BSgenome)){
      stop("If you state that custom_genome = T, you need to provide custom_location_RefCDS AND custom_BSgenome")
    }

  }


  MAF <- MAF_file
  if(length(which(MAF[,pos_column]==150713902))>0){
    MAF <- MAF[-which(MAF[,pos_column]==150713902),]
  }
  if(length(which(MAF[,pos_column]==41123095))>0){
    MAF <- MAF[-which(MAF[,pos_column]==41123095),]
  }
  #


  MAF <- MAF[,c(sample_ID_column,chr_column,pos_column,ref_column,alt_column)]

  message("Checking if any reference alleles provided do not match reference genome...")

  MAF <- MAF[,c(sample_ID_column,chr_column,pos_column,ref_column,alt_column)]



  if(custom_genome){
    reference_alleles <- as.character(BSgenome::getSeq(custom_BSgenome, paste("chr",MAF[,chr_column],sep=""),
                                                       strand="+", start=MAF[,pos_column], end=MAF[,pos_column]))

    if(length(which(MAF[,"Reference_Allele"] != reference_alleles))>0){
      message(paste(length(which(MAF[,"Reference_Allele"] != reference_alleles))," positions do not match out of ", nrow(MAF), ", (",round(length(which(MAF[,"Reference_Allele"] != reference_alleles))*100/nrow(MAF),4),"%), removing them from the analysis",sep=""))

      MAF <- MAF[-which(MAF[,"Reference_Allele"] != reference_alleles),]
    }

  }else{
    reference_alleles <- as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, paste("chr",MAF[,chr_column],sep=""),
                                                       strand="+", start=MAF[,pos_column], end=MAF[,pos_column]))

    if(length(which(MAF[,"Reference_Allele"] != reference_alleles))>0){
      message(paste(length(which(MAF[,"Reference_Allele"] != reference_alleles))," positions do not match out of ", nrow(MAF), ", (",round(length(which(MAF[,"Reference_Allele"] != reference_alleles))*100/nrow(MAF),4),"%), removing them from the analysis",sep=""))

      MAF <- MAF[-which(MAF[,"Reference_Allele"] != reference_alleles),]
    }
  }

  # message("Checking if any reference alleles provided do not match reference genome...")
  # MAF <- cancereffectsizeR::wrongRef_dndscv_checker(mutations = MAF)

  trinuc_data <- cancereffectsizeR::trinuc_profile_function_with_weights(
    input.MAF = MAF,
    sample.ID.column = sample_ID_column,
    chr.column = chr_column,
    pos.column = pos_column,
    ref.column = ref_column,
    alt.column = alt_column,
    bsg_for_deconstructSigs = custom_BSgenome)

  # this_cov_pca <- get(load(covariate_file))
  if(!is.null(covariate_file)){
    if(covariate_file=="dndscv_default"){
      data("covariates_hg19",package = "dndscv")
      genes_in_pca <- rownames(covs)
    }else{
      this_cov_pca <- get(data(list=covariate_file, package="cancereffectsizeR"))
      genes_in_pca <- rownames(this_cov_pca$rotation)
    }
  }

  path_to_library <- dir(.libPaths(),full.names=T)[grep(dir(.libPaths(),full.names=T),pattern="cancereffectsizeR")][1] # find the path to this package


  data("RefCDS_TP53splice",package = "cancereffectsizeR")

  if(custom_genome){
    dndscvout <- dndscv::dndscv(
      mutations = MAF,
      gene_list = genes_in_pca,
      cv = if(is.null(covariate_file)){NULL}else{if(covariate_file=="dndscv_default"){"hg19"}else{ this_cov_pca$rotation}},
      refdb = custom_location_RefCDS)
    load(custom_location_RefCDS)
  }else{
    dndscvout <- dndscv::dndscv(
      mutations = MAF,
      gene_list = genes_in_pca,
      cv = if(is.null(covariate_file)){NULL}else{if(covariate_file=="dndscv_default"){"hg19"}else{ this_cov_pca$rotation}},
      refdb = paste(path_to_library,"/data/RefCDS_TP53splice.RData",sep=""))
    data("refcds_hg19", package="dndscv")
  }




  RefCDS.our.genes <- RefCDS[which(sapply(RefCDS, function(x) x$gene_name) %in% dndscvout$genemuts$gene_name)]


  if(dndscvout$nbreg$theta>1){
    mutrates <- ((dndscvout$genemuts$n_syn + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv))/sapply(RefCDS.our.genes, function(x) colSums(x$L)[1]))/length(unique(MAF[,sample_ID_column]))
  }else{
    mutrates <- rep(NA,length(dndscvout$genemuts$exp_syn_cv))

    syn.sites <- sapply(RefCDS.our.genes, function(x) colSums(x$L)[1])

    for(i in 1:length(mutrates)){
      if( dndscvout$genemuts$exp_syn_cv[i] >  ((dndscvout$genemuts$n_syn[i] + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv[i])))){
        mutrates[i] <- (dndscvout$genemuts$exp_syn_cv[i]/syn.sites[i])/length(unique(MAF[,sample_ID_column]))
      }else{
        mutrates[i] <- (((dndscvout$genemuts$n_syn[i] + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv[i])))/syn.sites[i])/length(unique(MAF[,sample_ID_column]))
      }

    }
  }
  names(mutrates) <- dndscvout$genemuts$gene_name



  dndscv.pq <- dndscvout$sel_cv
  dndscv.pq$gene <- dndscv.pq$gene_name
  dndscv.pq$p <- dndscv.pq$pallsubs_cv
  dndscv.pq$q <- dndscv.pq$qallsubs_cv



  MAF_ranges <- GenomicRanges::GRanges(seqnames = MAF[,chr_column], ranges = IRanges::IRanges(start=MAF[,pos_column],end = MAF[,pos_column]))


  gene_name_matches <- GenomicRanges::nearest(x = MAF_ranges, subject = gr_genes,select=c("all"))

  MAF$Gene_name <- NA

  # take care of single hits
  single.choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))==1)])

  MAF$Gene_name[single.choice] <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)[which(S4Vectors::queryHits(gene_name_matches) %in% single.choice)]]


  multi.choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))>1)])

  # MAF[which(MAF$Start_Position==55118815),]

  all.possible.names <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)]
  query.spots <- S4Vectors::queryHits(gene_name_matches)

  for(i in 1:length(multi.choice)){
    # genes.for.this.choice <- gr_genes$names[subjectHits(gene_name_matches[which(queryHits(gene_name_matches)==multi.choice[i])])]
    genes.for.this.choice <- all.possible.names[which(query.spots==multi.choice[i])]
    if(length(which( genes.for.this.choice %in% names(mutrates) ))>0){
      MAF$Gene_name[multi.choice[i]] <- genes.for.this.choice[which( genes.for.this.choice %in% names(mutrates) )[1]]
    }else{
      MAF$Gene_name[multi.choice[i]] <- "Indeterminate"
    }
  }

  message("Calculating selection intensity...")



  selection_output <- cancereffectsizeR::selection_intensity_calculation(
    MAF_for_analysis = MAF,
    genes_for_analysis = genes_for_effect_size_analysis,
    mut_rates = mutrates,
    trinuc.mutation_data = trinuc_data$trinuc.mutation_data,
    dndscv_siggenes=dndscv.pq,
    output_from_mainMAF = output_from_mainMAF)

  return(list(selection_output=selection_output,
              mutation_rates=mutrates,
              trinuc_data=trinuc_data,dndscvout=dndscvout))

}



