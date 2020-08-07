#' Skip mutational signature analysis and assign group average relative trinucleotide-context-specific mutation rates to all samples
#'
#' This function calculates the relative rates of trinucleotide-context-specific mutations across
#' all SNV records from whole-exome and whole-genome MAF data and naively assigns these rates to all samples. 
#' This can be helpful if you do not have SNV mutational signatures available for your species, or if 
#' you want to assume that all samples share the same SNV mutational processes without relying on signatures.
#' Normally, if mutational signatures are available, it is better to use trinuc_snv_mutation_rates().
#' 
#' To reduce the influence of selection, only non-recurrent mutations (i.e., mutations that occur
#' in just one sample) are used to calculate the rates. Targeted sequencing data is excluded for
#' the same reason, and also because the trinucleotide composition of targeted regions could be
#' very different from that of the exome/genome.
#' 
#' @param cesa CESAnalysis object
#' @export
assign_group_average_trinuc_rates = function(cesa) {
  if(is.null(cesa) || ! is(cesa, "CESAnalysis")) {
    stop("Expected cesa to be a CESAnalysis object")
  }
  if(cesa@maf[, .N] == 0) {
    stop("No MAF data in the CESAnalysis")
  }
  if(all(cesa@samples$coverage == "targeted")) {
    stop("We can't estimate relative trinucleotide mutation rates without some exome/genome data in the CESAnalysis (all data is targeted sequencing).")
  }
  
  # Take just SNVs
  snv_maf = cesa@maf[Variant_Type == "SNV"]
  
  # Remove all recurrent SNVs (SNVs appearing in more than one sample)
  duplicated_vec_first <- duplicated(snv_maf[,.(Chromosome, Start_Position, Tumor_Allele)])
  duplicated_vec_last <- duplicated(snv_maf[,.(Chromosome, Start_Position, Tumor_Allele)],fromLast=T)
  duplicated_vec_pos <- which(duplicated_vec_first | duplicated_vec_last)
  if (length(duplicated_vec_pos) > 0) {
    snv_maf <- snv_maf[-duplicated_vec_pos,]
  }
  
  # Subset to just WES/WGS data
  non_tgs_samples = cesa@samples[coverage != "targeted", Unique_Patient_Identifier]
  snv_maf = snv_maf[Unique_Patient_Identifier %in% non_tgs_samples]
  
  
  # get trinuc contexts of each SNV and produce a data frame of counts, organized the same way as deconstructSigs data
  bsg = .ces_ref_data[[cesa@ref_key]]$genome
  trinuc = BSgenome::getSeq(bsg, snv_maf$Chromosome, snv_maf$Start_Position - 1, snv_maf$Start_Position + 1, as.character = T)
  
  # internal dict converts trinuc/mut (e.g., GTA:C) into deconstructSigs format ("G[T>C]A")
  ds_muts = factor(deconstructSigs_notations[.(trinuc, snv_maf$Tumor_Allele), deconstructSigs_ID], levels = deconstructSigs_trinuc_string)
  
  # mysteriously convert two-way table to data frame
  tmp = table(snv_maf$Unique_Patient_Identifier, ds_muts)
  counts = apply(tmp, 2, rbind)
  rownames(counts) = rownames(tmp)
  trinuc_breakdown_per_tumor = as.data.frame(counts)
  
  # produce normalized rates (putting in pseudocounts any are 0)
  trinuc_prop = colSums(trinuc_breakdown_per_tumor) / sum(trinuc_breakdown_per_tumor) 
  
  if(any(trinuc_prop == 0)) {
    lowest_nonzero_rate = min(trinuc_prop[trinuc_prop != 0])
    trinuc_prop = trinuc_prop + lowest_nonzero_rate
    trinuc_prop = trinuc_prop / sum(trinuc_prop) # renormalizing
  }
  
  # create trinuc_proportion_matrix (1 row per sample, all rows with identical trinuc_prop)
  num_samples = cesa@samples[, .N]
  trinuc_proportion_matrix = matrix(rep(trinuc_prop, num_samples), byrow = T, ncol = 96, 
                                    dimnames = list(cesa@samples$Unique_Patient_Identifier, names(trinuc_prop)))

  cesa@status[["trinucleotide mutation rates"]] = "from assign_group_average_trinuc_rates"
  cesa@status[["gene mutation rates"]] = "uncalculated (run gene_mutation_rates)"
  cesa@trinucleotide_mutation_weights = list(trinuc_proportion_matrix=trinuc_proportion_matrix)
  return(cesa)
}