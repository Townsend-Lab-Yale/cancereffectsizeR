#' \code{\link[dndscv]{dndscv}} wrong reference checker
#'
#' Checks if reference base does not match reference within \code{\link[dndscv]{dndscv}} and tells you if there is a mismatch so you can remove. All inputs and most code same as in \code{\link[dndscv]{dndscv}}. See ?dndscv for more information.
#'
#' @param mutations Same as in \code{\link[dndscv]{dndscv}}. Table of
#'  mutations (5 columns: sampleID, chr, pos, ref, alt). Only list independent events as mutations.
#' @return matrix of provided mutations that do not match reference genome
#' @import dndscv
#' @import seqinr
#' @import Biostrings
#' @import MASS
#' @import GenomicRanges
#' @export


wrongRef_dndscv_checker = function(mutations, gene_list = NULL, refdb = "hg19", sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = 3, max_coding_muts_per_sample = 3000, use_indel_sites = T, min_indels = 5, maxcovs = 20, constrain_wnon_wspl = T, outp = 3) {

  original_mutation_input <- mutations
  ## 1. Environment
  message("[1] Loading the environment...")

  mutations[,c(1,2,4,5)] = lapply(mutations[,c(1,2,4,5)], as.character) # Factors to character

  # [Input] Reference database
  if (refdb == "hg19") {
    data("refcds_hg19", package="dndscv")
  } else {
    load(refdb)
  }

  # [Input] Gene list (The user can input a gene list as a character vector)
  if (is.null(gene_list)) {
    gene_list = sapply(RefCDS, function(x) x$gene_name) # All genes [default]
  } else { # Using only genes in the input gene list
    allg = sapply(RefCDS,function(x) x$gene_name)
    nonex = gene_list[!(gene_list %in% allg)]
    if (length(nonex)>0) { stop(sprintf("The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", "))) }
    RefCDS = RefCDS[allg %in% gene_list] # Only input genes
    gr_genes = gr_genes[gr_genes$names %in% gene_list] # Only input genes
  }

  # [Input] Covariates (The user can input a custom set of covariates as a matrix)
  if (is.character(cv)) {
    data(list=sprintf("covariates_%s",cv), package="dndscv")
  } else {
    covs = cv
  }

  # [Input] Known cancer genes (The user can input a gene list as a character vector)
  if (kc[1] %in% c("cgc81")) {
    data(list=sprintf("cancergenes_%s",kc), package="dndscv")
  } else {
    known_cancergenes = kc
  }

  # [Input] Substitution model (The user can also input a custom substitution model as a matrix)
  if (length(sm)==1) {
    data(list=sprintf("submod_%s",sm), package="dndscv")
  } else {
    substmodel = sm
  }

  # Expanding the reference sequences [for faster access]
  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
    RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
    RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
    if (!is.null(RefCDS[[j]]$seq_splice)) {
      RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
      RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
      RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
    }
  }


  ## 2. Mutation annotation
  message("[2] Annotating the mutations...")

  colnames(mutations) = c("sampleID","chr","pos","ref","mut")
  nt = c("A","C","G","T")
  trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  trinucinds = setNames(1:64, trinucs)

  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)

  ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
  gr_genes_ind = ind[gr_genes$names]

  # Warning about possible unannotated dinucleotide substitutions
  if (any(diff(mutations$pos)==1)) {
    warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
  }

  # Warning about multiple instances of the same mutation in different sampleIDs
  if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
    warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
  }

  # Mapping mutations to genes
  gr_muts = GenomicRanges::GRanges(mutations$chr, IRanges::IRanges(mutations$pos,mutations$pos))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene
  mutations$geneind = gr_genes_ind[ol[,2]]
  mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]

  # Optional: Excluding samples exceeding the limit of mutations/sample [see Default parameters]
  nsampl = sort(table(mutations$sampleID))
  exclsamples = NULL
  if (any(nsampl>max_coding_muts_per_sample)) {
    message(sprintf('    Note: %0.0f samples excluded for exceeding the limit of mutations per sample',sum(nsampl>max_coding_muts_per_sample)))
    exclsamples = names(nsampl[nsampl>max_coding_muts_per_sample])
    mutations = mutations[!(mutations$sampleID %in% names(nsampl[nsampl>max_coding_muts_per_sample])),]
  }

  # Optional: Limiting the number of mutations per gene per sample (to minimise the impact of unannotated kataegis and other mutation clusters) [see Default parameters]
  mutrank = ave(mutations$pos, paste(mutations$sampleID,mutations$gene), FUN = function(x) rank(x))
  exclmuts = NULL
  if (any(mutrank>max_muts_per_gene_per_sample)) {
    message(sprintf('    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample',sum(mutrank>max_muts_per_gene_per_sample)))
    exclmuts = mutations[mutrank>max_muts_per_gene_per_sample,]
    mutations = mutations[mutrank<=max_muts_per_gene_per_sample,]
  }

  # Additional annotation of substitutions

  snv = (mutations$ref %in% nt & mutations$mut %in% nt)
  indels = mutations[!snv,]
  mutations = mutations[snv,]
  mutations$ref_cod = mutations$ref
  mutations$mut_cod = mutations$mut
  compnt = setNames(rev(nt), nt)

  mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
  isminus = (mutations$strand==-1)
  mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
  mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]

  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
  }

  # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals

  chr2cds = function(pos,cds_int,strand) {
    if (strand==1) {
      return(which(pos==unlist(apply(cds_int, 1, function(x) x[1]:x[2]))))
    } else if (strand==-1) {
      return(which(pos==rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2])))))
    }
  }

  # Annotating the functional impact of each substitution and populating the N matrices

  ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = array(NA, nrow(mutations))

  for (j in 1:nrow(mutations)) {

    geneind = mutations$geneind[j]
    pos = mutations$pos[j]
    if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution

      impact[j] = "Essential_Splice"; impind = 4
      pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
      cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      aachange[j] = ntchange[j] = "-"

    } else { # Coding substitution

      pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
      cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
      old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
      pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
      new_codon = old_codon; new_codon[pos_in_codon] = mutations$mut_cod[j]
      old_aa = seqinr::translate(old_codon)
      new_aa = seqinr::translate(new_codon)
      aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
      ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])

      # Annotating the impact of the mutation
      if (new_aa == old_aa){
        impact[j] = "Synonymous"; impind = 1
      } else if (new_aa == "*"){
        impact[j] = "Nonsense"; impind = 3
      } else if (old_aa != "*"){
        impact[j] = "Missense"; impind = 2
      } else if (old_aa=="*") {
        impact[j] = "Stop_loss"; impind = NA
      }
    }

    if (mutations$ref_cod[j] != as.character(cdsnt)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
      wrong_ref[j] = 1
    } else if (!is.na(impind)) { # Correct base annotation in the input mutation file
      trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
      RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices
    }

    if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g %%...', round(j/nrow(mutations),2)*100)) }
  }

  mutations$ref3_cod = ref3_cod
  mutations$mut3_cod = mut3_cod
  mutations$aachange = aachange
  mutations$ntchange = ntchange
  mutations$impact = impact
  mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]

  # return(mutations[which(wrong_ref==1),])

  if(nrow(mutations[which(wrong_ref==1),])>0){
    message(paste(nrow(wrong_refs),"mismatched mutations, removing them and returning remaining mutations."))
    for(i in 1:nrow(wrong_refs)){
      if(length(which(original_mutation_input$Chromosome==wrong_refs$chr[i] & original_mutation_input$Start_Position==wrong_refs$pos[i] &
                      original_mutation_input$Unique_patient_identifier==wrong_refs$sampleID[i] & original_mutation_input$Tumor_allele==wrong_refs$mut[i]))>0){
        original_mutation_input <- original_mutation_input[-which(original_mutation_input$Chromosome==wrong_refs$chr[i] & original_mutation_input$Start_Position==wrong_refs$pos[i] &
                            original_mutation_input$Unique_patient_identifier==wrong_refs$sampleID[i] & original_mutation_input$Tumor_allele==wrong_refs$mut[i]),]
      }
    }
    return(original_mutation_input)
  }else{
    message("No mismatched mutations! Returning original input.")
   return(original_mutation_input)
  }

} # EOF
