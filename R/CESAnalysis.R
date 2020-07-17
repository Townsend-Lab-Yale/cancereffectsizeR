#' Create a cancereffectsizeR analysis
#' 
#' @description Creates a CESAnalysis object, the central data structure of cancereffectsizeR
#' @param progression_order evolutionary order of tumor stage progression (e.g. c("Primary", "Metastatic"))
#' @param genome Genome build of MAF data (currently, just hg19 supported)
#' @return CESAnalysis object
#' @export
CESAnalysis = function(genome = NULL, progression_order = NULL) {
  genome_dir = get_genome_data_directory(genome)
  genome_path = paste0(genome_dir, "/genome_package_name.rds")
  if(! file.exists(genome_path)) {
    stop(paste0("Something is wrong with the genome data installation.\n",
                "Expected to find a reference genome at ", genome_path, "."))
  }
  genome_package = readRDS(genome_path)
  genome = BSgenome::getBSgenome(genome_package)
  GenomeInfoDb::seqlevelsStyle(genome) = "NCBI"
  message(crayon::black(paste0("Okay, this CES analysis will use the ", 
                               tolower(BSgenome::commonName(genome)), " genome (", BSgenome::releaseName(genome), ").")))
  message(crayon::black("Note: We'll be using NCBI-style chromosome names (i.e., no \"chr\" prefixes)."))
  
  
  # Validate progression_order
  if (is.null(progression_order)) {
    progression_order = c("stageless")
  }
  
  if (is.numeric(progression_order)) {
    progression_order = as.character(progression_order)
  }

  if (! is(progression_order, "character")) {
    stop("progression_order should be a character vector of chronological tumor states (e.g., Primary, Metastatic)")
  }

  
  # simple analysis status tracking; used to guide user in show(CESAnalysis)
  status = list("genome" = GenomeInfoDb::providerVersion(genome),
                "progressions" = paste0(progression_order, collapse = ", "))
  if(length(progression_order) == 1) {
    status[["progressions"]] = NULL
  }
  advanced = list("version" = packageVersion("cancereffectsizeR"))
  cesa = new("CESAnalysis", status = status, genome = genome, maf = data.table(), excluded = data.table(),
             progressions = progression_order, mutrates = data.table(),
             gene_epistasis_results = data.table(), selection_results = data.table(), genome_data_dir = genome_dir,
             advanced = advanced, samples = data.table(), mutations = list())
  return(cesa)
}

#' View data loaded into CESAnalysis
#' 
#' returns a data.table containing MAF records used in the given CESAnalysis
#' @param cesa CESAnalysis object
#' @export
maf_records = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: maf_records(cesa), where cesa is a CESAnalysis")
  }
  if(cesa@maf[,.N] == 0) {
    stop("No MAF data has been loaded")
  }
  return(cesa@maf)
}

#' View excluded MAF data
#' 
#' returns a data.table containing MAF records that were excluded from the given CESAnalysis
#' @param cesa CESAnalysis object
#' @export
excluded_maf_records = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: excluded_maf_records(cesa), where cesa is a CESAnalysis")
  }
  if(cesa@maf[,.N] == 0) {
    stop("No MAF data has been loaded yet, so naturally no records have been excluded.")
  }
  if(cesa@excluded[,.N] == 0) {
    message("Returned an empty data table since no records have been excluded.")
  }
  return(cesa@excluded)
}

#' View sample metadata
#' 
#' returns a data.table with info on all samples in the CESAnalysis
#' @param cesa CESAnalysis object
#' @export
get_sample_info = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_sample_info(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@samples[,.N] == 0) {
    stop("No MAF data has been loaded yet, so naturally there is no sample data.")
  }
  
  # user doesn't need progression columns for single-progression-state analyses
  if(length(cesa@progressions) == 1) {
    return(cesa@samples[, -c("progression_index", "progression_name")])
  } else {
    return(cesa@samples)
  }
  return(cesa@samples)
}

#' Get expected relative trinucleotide-specific SNV mutation rates
#' 
#' @param cesa CESAnalysis object
#' @export
get_trinuc_rates = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_trinuc_rates(cesa), where cesa is a CESAnalysis")
  }
  return(as.data.table(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix, keep.rownames = "Unique_Patient_Identifier"))
}

#' Get table of signature weights by tumor
#' 
#' @param cesa CESAnalysis object
#' @param include_tumors_without_data include rows (consisting of group-average weights) for samples that did not undergo
#'                                    signature extraction, such as TGS tumors or the rare sample with no non-recurrent SNVs
#' @export
get_signature_weights = function(cesa = NULL, include_tumors_without_data = FALSE) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_signature_weights(cesa), where cesa is a CESAnalysis")
  }
  sig_table = cesa@trinucleotide_mutation_weights$signature_weight_table
  
  if (include_tumors_without_data) {
    tumors_without_data = setdiff(cesa@samples$Unique_Patient_Identifier, sig_table$Unique_Patient_Identifier)
    num_to_add = length(tumors_without_data)
    if (num_to_add > 0) {
      group_avg_weights = as.numeric(cesa@trinucleotide_mutation_weights$group_average_dS_output$adjusted_sig_output$weights)
      new_rows = matrix(nrow = num_to_add, data = rep.int(group_avg_weights, num_to_add), byrow = T)
      colnames(new_rows) = colnames(cesa@trinucleotide_mutation_weights$group_average_dS_output$adjusted_sig_output$weights)
      total_snvs = cesa@maf[Variant_Type == "SNV"][, .N, keyby = "Unique_Patient_Identifier"][tumors_without_data, N]
      total_snvs[is.na(total_snvs)] = 0
      new_table = data.table(Unique_Patient_Identifier = tumors_without_data, total_snvs = total_snvs, 
                             sig_extraction_snvs = 0, group_avg_blended = T)
      new_table = cbind(new_table, new_rows)
      sig_table = rbind(sig_table, new_table)
    }
  }
  return(sig_table)
}

#' Get table of neutral gene mutation rates by progression state
#' 
#' @param cesa CESAnalysis object
#' @export
get_gene_rates = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_gene_rates(cesa), where cesa is a CESAnalysis")
  }
  return(cesa@mutrates)
}

#' Get lists of mutations and annotations
#' 
#' @param cesa CESAnalysis object
#' @export
get_mutations = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: get_mutations(cesa), where cesa is a CESAnalysis")
  }
  return(cesa@mutations)
}



#' View results from ces_snv
#' 
#' returns a data table of SNV effect sizes generated with ces_snv
#' @param cesa CESAnalysis object
#' @export
snv_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: snv_results(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@selection_results[, .N] == 0) {
    stop("No results yet from ces_snv in this CESAnalysis")
  }
  return(cesa@selection_results)
}

#' View results from gene-level epistasis analysis
#' 
#' returns a data table of pairwise gene epistasis effect sizes generated with ces_gene_epistasis
#' @param cesa CESAnalysis object
#' @export
gene_epistasis_results = function(cesa = NULL) {
  if(! is(cesa, "CESAnalysis")) {
    stop("\nUsage: gene_epistasis_results(cesa), where cesa is a CESAnalysis")
  }
  if (cesa@gene_epistasis_results[, .N] == 0) {
    stop("No results yet from ces_gene_epistasis in this CESAnalysis")
  }
  return(cesa@gene_epistasis_results)
}

