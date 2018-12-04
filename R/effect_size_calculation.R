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
#' @param covariate_file Either NULL and uses the \code{dndscv}
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
                            genes_for_effect_size_analysis="all"){

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
  MAF <- cancereffectsizeR::wrongRef_dndscv_checker(mutations = MAF)

  trinuc_data <- cancereffectsizeR::trinuc_profile_function_with_weights(
    input_MAF = MAF,
    sample_ID_column = sample_ID_column,
    chr_column = chr_column,
    pos_column = pos_column,
    ref_column = ref_column,
    alt_column = alt_column)

  # this_cov_pca <- get(load(covariate_file))

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
    mutations = MAF,
    gene_list = genes_in_pca,
    cv = if(is.null(covariate_file)){ "hg19"}else{ this_cov_pca$rotation},
    refdb = paste(path_to_library,"/data/RefCDS_TP53splice.RData",sep=""))



  data("refcds_hg19", package="dndscv")

  RefCDS_our_genes <- RefCDS[which(sapply(RefCDS, function(x) x$gene_name) %in% dndscvout$genemuts$gene_name)]


  if(dndscvout$nbreg$theta>1){
    mutrates <- ((dndscvout$genemuts$n_syn + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv))/sapply(RefCDS_our_genes, function(x) colSums(x$L)[1]))/length(unique(MAF[,sample_ID_column]))
  }else{
    mutrates <- rep(NA,length(dndscvout$genemuts$exp_syn_cv))

    syn_sites <- sapply(RefCDS_our_genes, function(x) colSums(x$L)[1])

    for(i in 1:length(mutrates)){
      if( dndscvout$genemuts$exp_syn_cv[i] >  ((dndscvout$genemuts$n_syn[i] + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv[i])))){
        mutrates[i] <- (dndscvout$genemuts$exp_syn_cv[i]/syn_sites[i])/length(unique(MAF[,sample_ID_column]))
      }else{
        mutrates[i] <- (((dndscvout$genemuts$n_syn[i] + dndscvout$nbreg$theta - 1)/(1+(dndscvout$nbreg$theta/dndscvout$genemuts$exp_syn_cv[i])))/syn_sites[i])/length(unique(MAF[,sample_ID_column]))
      }

    }
  }
  names(mutrates) <- dndscvout$genemuts$gene_name



  dndscv_pq <- dndscvout$sel_cv
  dndscv_pq$gene <- dndscv_pq$gene_name
  dndscv_pq$p <- dndscv_pq$pallsubs_cv
  dndscv_pq$q <- dndscv_pq$qallsubs_cv



  MAF_ranges <- GenomicRanges::GRanges(seqnames = MAF[,chr_column], ranges = IRanges::IRanges(start=MAF[,pos_column],end = MAF[,pos_column]))


  gene_name_matches <- GenomicRanges::nearest(x = MAF_ranges, subject = gr_genes,select=c("all"))

  MAF$Gene_name <- NA

  # take care of single hits
  single_choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))==1)])

  MAF$Gene_name[single_choice] <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)[which(S4Vectors::queryHits(gene_name_matches) %in% single_choice)]]


  multi_choice <- as.numeric(names(table(S4Vectors::queryHits(gene_name_matches)))[which(table(S4Vectors::queryHits(gene_name_matches))>1)])

  # MAF[which(MAF$Start_Position==55118815),]

  all_possible_names <- gr_genes$names[S4Vectors::subjectHits(gene_name_matches)]
  query_spots <- S4Vectors::queryHits(gene_name_matches)

  for(i in 1:length(multi_choice)){
    # genes.for.this.choice <- gr_genes$names[subjectHits(gene_name_matches[which(queryHits(gene_name_matches)==multi.choice[i])])]
    genes_for_this_choice <- all_possible_names[which(query_spots==multi_choice[i])]
    if(length(which( genes_for_this_choice %in% names(mutrates) ))>0){
      MAF$Gene_name[multi_choice[i]] <- genes_for_this_choice[which( genes_for_this_choice %in% names(mutrates) )[1]]
    }else{
      MAF$Gene_name[multi_choice[i]] <- "Indeterminate"
    }
  }

  message("Calculating selection intensity...")



  selection_output <- cancereffectsizeR::selection_intensity_calculation(
    MAF_for_analysis = MAF,
    genes_for_analysis = genes_for_effect_size_analysis,
    mut_rates = mutrates,
    trinuc_mutation_data = trinuc_data$trinuc.mutation_data,
    dndscv_siggenes=dndscv_pq,
    output_from_mainMAF = output_from_mainMAF)

  return(list(selection_output=selection_output,
              mutation_rates=mutrates,
              trinuc_data=trinuc_data,dndscvout=dndscvout))

}



