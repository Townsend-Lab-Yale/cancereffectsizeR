#' Baseline mutation rate calculation
#' 
#' Caculates neutral mutation rates at specific sites based on gene mutation rates and the relative
#' trinucleotide-context-specific SNV mutation rates of each sample
#' 
#' @param cesa CESAnalysis with gene mutation rates and tumor-specific trinucleotide-context-specific mutation rates already calculated
#' @param aac_ids vector of IDs for amino acid change variants
#' @param snv_ids vector of IDs for SNVs
#' @param variant_ids vector of mixed IDs (faster to use snv_ids and aac_ids for large jobs, if already known)
#' @param samples vector of sample IDs (Unique_Patient_Identifier) to include in mutation rate table (defaults to all samples)
#' @param cores number of cores to use for mutation processing (useful for large data sets or mutation lists)
#' @return a data table of mutation rates with one column per variant, and a Unique_Patient_Identifier column identifying each row
#' @export
baseline_mutation_rates = function(cesa, aac_ids = NULL, snv_ids = NULL, variant_ids = NULL, samples = NULL, cores = 1) {
  
  if(! cesa@advanced$trinuc_done) {
    stop("Some samples lack trinucleotide-context-specific mutation rates, so site-level mutation rates can't be calculated yet.")
  }
  if(! cesa@advanced$gene_rates_done) {
    stop("Some samples lack gene mutation rates, so site-level mutation rates can't be calculated yet.")
  }
  
  if(is.null(cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix)) {
    stop("No trinucleotide mutation rates found, so can't calculate variant-level mutation rates.")
  }
  
  if (! is.null(variant_ids) && (! is.null(snv_ids) || ! is.null(aac_ids))) {
    stop("You can use snv_ids/aac_ids or variant_ids, but not a combination.")
  }
  
  if (! is.null(aac_ids)) {
    if(! is.character(aac_ids)) {
      stop("aac_ids should be type character")
    }
  }
  
  if (! is.null(snv_ids)) {
    if(! is.character(snv_ids)) {
      stop("snv_ids should be type character")
    }
  }
  
  if (! is.null(variant_ids)) {
    if(! is.character(variant_ids)) {
      stop("variant_ids should be type character")
    }
  }
  
  # Helping the user out
  aac_ids = unique(na.omit(aac_ids))
  snv_ids = unique(na.omit(snv_ids))
  variant_ids = unique(na.omit(variant_ids))
  
  # If variant IDs is used, snv/aac IDs were not
  if (! is.null(variant_ids)) {
    aac_ids = cesa@mutations$amino_acid_change[variant_ids, aac_id, nomatch = NULL]
    snv_ids = cesa@mutations$snv[variant_ids, snv_id, nomatch = NULL]
    
    if (length(snv_ids) + length(aac_ids) != length(variant_ids)) {
      unknown_ids = setdiff(variant_ids, union(aac_ids, snv_ids))
      stop(length(unknown_ids), " variant IDs are either invalid or not present in annotation tables: ", 
           paste(unknown_ids, collapse = ", "), ".")
    }
  } else {
    missing_aac = setdiff(aac_ids, cesa@mutations$amino_acid_change[aac_ids, aac_id, nomatch = NULL])
    missing_snv = setdiff(snv_ids, cesa@mutations$snv[snv_ids, snv_id, nomatch = NULL])
    all_missing = c(missing_snv, missing_aac)
    num_missing = length(all_missing)
    if(num_missing > 0) {
      stop(num_missing, " variant IDs are either invalid or not present in annotation tables: ", 
           paste(all_missing, collapse = ", "), ".")
    }
  }
  
  # Let user specify a subset of samples to calculate rates (or, by default, use all samples)
  if(! is.null(samples)) {
    if (! is.character(samples)) {
      stop("Samples should be character vector")
    }
    samples = unique(na.omit(samples))
    subsetted_samples = cesa@samples[samples, on = "Unique_Patient_Identifier", nomatch = NULL]
    if(subsetted_samples[, .N] != length(samples)) {
      invalid_samples = setdiff(samples, subsetted_samples$Unique_Patient_Identifier)
      stop(length(samples), " samples are not in the CESAnalysis samples table: ", 
           paste(invalid_samples, collapse = ", "), ".")
    }
    samples = subsetted_samples
  } else {
    samples = cesa@samples
  }
  mutations = cesa@mutations
  
  # can drop AAC mutations not requested
  if(length(aac_ids) > 0) {
    mutations$amino_acid_change = mutations$amino_acid_change[aac_ids]
    setkey(mutations$amino_acid_change, "aac_id")
  } else {
    mutations$amino_acid_change = data.table()
  }
  
  
  # Give a progress message if this is going to take more than a few seconds
  num_variants = length(aac_ids) + length(snv_ids)
  
  if(num_variants == 0) {
    stop("Can't calculate mutation rates because no variants were input.")
  }
  
  num_samples = samples[, .N]
  if(num_variants * num_samples > 1e7) {
    num_variants = format(num_variants, big.mark = ",") # big.mark? R is so weird
    num_samples = format(num_samples, big.mark = ",")
    message(sprintf("Preparing to calculate baseline mutation rates in %s samples across %s sites...", num_samples, num_variants))
  }
  
  # produce a table with all pairwise combinations of Unique_Patient_Identifier and relevant regional rates
  # relevant genes/pids are those associated with one of the AACs/SNVs of interest
  using_pid = FALSE
  if ("pid" %in% names(cesa@mutrates)) {
    using_pid = TRUE
    # Note that unlist(NULL data.table) returns NULL (which is okay)
    relevant_regions = union(mutations$amino_acid_change$pid, unlist(mutations$snv[snv_ids, nearest_pid, on = "snv_id"]))
  } else {
    relevant_regions = union(mutations$amino_acid_change$gene, unlist(mutations$snv[snv_ids, genes, on = "snv_id"]))
  }
  sample_region_rates = as.data.table(expand.grid(region = relevant_regions, Unique_Patient_Identifier = samples$Unique_Patient_Identifier, 
                                                  stringsAsFactors = F), key = "Unique_Patient_Identifier")

  # add gene (regional) mutation rates to the table by using @mutrates and the groups of each samples
  sample_region_rates = sample_region_rates[samples[, .(Unique_Patient_Identifier, gene_rate_grp)]]
  
  if (using_pid) {
    # there's also a gene given in transcript rates tables, but we'll drop it here
    melted_mutrates = melt.data.table(cesa@mutrates[relevant_regions, -"gene", on = "pid"], id.vars = c("pid"))
    setnames(melted_mutrates, "pid", "region")
  } else {
    melted_mutrates = melt.data.table(cesa@mutrates[relevant_regions, on = "gene"], id.vars = c("gene"))
    setnames(melted_mutrates, "gene", "region")
  }
  
  setnames(melted_mutrates, c("variable", "value"), c("gene_rate_grp", "raw_rate"))
  sample_region_rates = melted_mutrates[sample_region_rates, , on = c("region", "gene_rate_grp")]
  
  cds_trinuc_comp = get_ref_data(cesa, "gene_trinuc_comp")
  
  
  # Hash trinuc rates for faster runtime with large data sets (where there could be millions of queries of trinuc_mat)
  # Will also be using the trinuc_mat directly for summed rates, so need to subset it to just samples of interest
  trinuc_rates = new.env(parent = emptyenv())
  trinuc_mat = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix[samples$Unique_Patient_Identifier, , drop = F]
  for (row in rownames(trinuc_mat)) { 
    trinuc_rates[[row]] = unname(trinuc_mat[row, ])
  }
  # dot product of trinuc comp and patient's expected relative trinuc rates yields the denominator
  # for site-specific mutation rate calculation; numerator is raw gene rate multipled by the patient's relative rate for the site's trinuc context
  # (this last value gets multipled in by get_baseline functions below)
  sample_region_rates[, aggregate_rate := raw_rate / sum(cds_trinuc_comp[[region]] * trinuc_rates[[Unique_Patient_Identifier]]), by = c("region", "Unique_Patient_Identifier")]
  setkey(sample_region_rates, "region")
  
  if(length(aac_ids > 0)) {
    trinuc_mut_by_aac = mutations$amino_acid_change[, .(trinuc_mut = list(mutations$snv[unlist(constituent_snvs), trinuc_mut])), by = "aac_id"]
    
    if (using_pid) {
      mutations$amino_acid_change[, region := pid]
    } else {
      mutations$amino_acid_change[, region := gene]
    }
    region_by_aac = mutations$amino_acid_change[aac_ids, region]
    aac_regions = unique(region_by_aac)
    
    # build vector agg_rates (see above comment) corresponding to each aac
    tmp = sample_region_rates[aac_regions, .(list(aggregate_rate)), by = c("region")]
    agg_rates_by_region = tmp$V1
    names(agg_rates_by_region) = tmp$region
    agg_rates_by_region = agg_rates_by_region[region_by_aac]
    agg_rates_by_region = lapply(agg_rates_by_region, unlist)
    
    # for each sample, multiply agg_rate by relative trinuc rate of the context for all AAC sites
    get_baseline_aac = function(aac, agg_rates) {
      trinuc_mut = trinuc_mut_by_aac[aac, unlist(trinuc_mut)]
      sample_rates = rowSums(trinuc_mat[, trinuc_mut, drop = F])
      return(as.numeric(agg_rates * sample_rates))
    }
    if(length(aac_ids) > 1000) {
      message("Calculating baseline rates for amino-acid-changing mutations...")
      pbopts = pbapply::pboptions()
    } else {
      pbopts = pbapply::pboptions(type = "none")
    }
    
    # data.table is supposed to automatically go to single-thread mode when running parallel lapply,
    # but apparently doesn't work in pbapply
    if (cores > 1) {
      original_dt_threads = setDTthreads(1)
    }
    aac_rate_list = pbapply::pblapply(1:length(aac_ids), function(x) get_baseline_aac(aac_ids[x], agg_rates_by_region[[x]]), cl = cores)
    if (cores > 1) {
      setDTthreads(original_dt_threads)
    }
    pbapply::pboptions(pbopts)
    names(aac_rate_list) = aac_ids
  } else {
    aac_rate_list = NULL
  }
  

  # repeat with SNVs (slightly different handling since SNVs can have more than one gene/pid associated)
  if (length(snv_ids) > 0) {
    trinuc_mut_by_snv = mutations$snv[snv_ids, trinuc_mut, keyby = "snv_id"]
    if (using_pid) {
      region_by_snv = mutations$snv[snv_ids, .(region = nearest_pid), by = "snv_id"]
    } else {
      region_by_snv = mutations$snv[snv_ids, .(region = genes), by = "snv_id"]
    }
    snv_regions = unique(unlist(region_by_snv$region))
    tmp = sample_region_rates[snv_regions, .(list(aggregate_rate)), by = "region"]
    agg_rates_by_region = tmp$V1
    names(agg_rates_by_region) = tmp$region
    
    get_baseline_snv = function(snv, region) {
      trinuc_mut = trinuc_mut_by_snv[snv, trinuc_mut]
      sample_rates = trinuc_mat[, trinuc_mut, drop = F]
      if(length(region) > 1) {
        agg_rates = rowMeans(simplify2array(agg_rates_by_region[region]))
      } else {
        agg_rates = agg_rates_by_region[[region]]
      }
      return(as.numeric(agg_rates * sample_rates))
    }
    if(length(snv_ids) > 1000) {
      message("Calculating baseline rates for noncoding SNVs...")
      pbopts = pbapply::pboptions()
    } else {
      # turn of progress bar when few variants
      pbopts = pbapply::pboptions(type = "none")
    }
    
    if (cores > 1) {
      original_dt_threads = setDTthreads(1)
    }
    snv_rate_list = pbapply::pblapply(1:length(snv_ids), function(x) get_baseline_snv(snv_ids[x], region_by_snv$region[[x]]), cl = cores)
    if (cores > 1) {
      setDTthreads(original_dt_threads)
    }
    pbapply::pboptions(pbopts)
    names(snv_rate_list) = snv_ids
  } else {
    snv_rate_list = NULL
  }
  
  
  # Combine all rates and create table
  baseline_rates = as.data.table(c(aac_rate_list, snv_rate_list))
  baseline_rates$Unique_Patient_Identifier = samples$Unique_Patient_Identifier
  setcolorder(baseline_rates, "Unique_Patient_Identifier")
  return(baseline_rates)
}


