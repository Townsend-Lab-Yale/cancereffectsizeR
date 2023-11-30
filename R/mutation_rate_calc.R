#' Baseline mutation rate calculation
#' 
#' Calculates neutral mutation rates at specific sites based on gene mutation rates and the relative
#' trinucleotide-context-specific SNV mutation rates of each sample
#' 
#' @param cesa CESAnalysis with gene mutation rates and tumor-specific trinucleotide-context-specific mutation rates already calculated
#' @param aac_ids vector of IDs for amino acid change variants
#' @param snv_ids vector of IDs for SNVs
#' @param variant_ids vector of mixed IDs (faster to use snv_ids and aac_ids for large jobs, if already known)
#' @param samples Which samples to calculate rates for. Defaults to all samples. Can be a
#'   vector of Unique_Patient_Identifiers, or a data.table containing rows from the
#'   CESAnalysis sample table.
#' @return a data table of mutation rates with one column per variant, and a Unique_Patient_Identifier column identifying each row
#' @export
baseline_mutation_rates = function(cesa, aac_ids = NULL, snv_ids = NULL, variant_ids = NULL, samples = character()) {
  
  if(! is(cesa, "CESAnalysis")) {
    stop("cesa should be CESAnalysis.")
  }

  samples = select_samples(cesa, samples)
  if(samples[, anyNA(sig_analysis_grp)]) {
    stop("Some input samples lack trinucleotide-context-specific relative mutation rates, so site-specific mutation rates can't be calculated yet.")
  }
  if(samples[, anyNA(gene_rate_grp)]) {
    stop("Some input samples lack gene mutation rates, so site-level mutation rates can't be calculated yet.")
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
  mutations = copy(cesa@mutations)
  
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
    num_variants = format(num_variants, big.mark = ",")
    num_samples = format(num_samples, big.mark = ",")
    message(sprintf("Preparing to calculate baseline mutation rates in %s samples across %s sites...", num_samples, num_variants))
  }
  # Get a copy of regional mutation rates, leaving out CI columns
  mutrates = cesa@mutrates[, .SD, .SDcols = patterns('gene|pid|(rate_grp_\\d+)$')]
  
  
  # produce a table with all pairwise combinations of Unique_Patient_Identifier and relevant regional rates
  # relevant genes/pids are those associated with one of the AACs/SNVs of interest
  using_pid = cesa@advanced$cds_refset
  if(using_pid) {
    if (! "pid" %in% names(mutrates)) {
      # If using a CDS-based refset (i.e., not ces.refset.hg19) and regional rates are gene-based
      # (this may happen if user inputs custom rates), add protein annotations to mutrates table.
      # All proteins associated with a given gene will get the same rate.
      gr_genes = get_ref_data(cesa, "gr_genes")
      if("gene" %in% names(GenomicRanges::mcols(gr_genes))) {
        pid_to_gene = unique(as.data.table(gr_genes)[, .(pid = names, gene = gene)])
      }
      mutrates = mutrates[pid_to_gene, on = 'gene']
    }
    # Note that unlist(NULL data.table) returns NULL (which is okay)
    relevant_regions = union(mutations$amino_acid_change$pid, unlist(mutations$snv[snv_ids, nearest_pid, on = "snv_id"]))
  } else {
    relevant_regions = union(mutations$amino_acid_change$gene, unlist(mutations$snv[snv_ids, genes, on = "snv_id"]))
  }
  
  # Get all combinations of genes (regions) and patients, then merge in gene_rate_grp from samples table.
  sample_region_rates = as.data.table(expand.grid(region = relevant_regions, Unique_Patient_Identifier = samples$Unique_Patient_Identifier, 
                                                  stringsAsFactors = F), key = "Unique_Patient_Identifier")
  sample_region_rates[samples, gene_rate_grp := paste0('rate_grp_', gene_rate_grp), on = 'Unique_Patient_Identifier']
  
  if (using_pid) {
    # there's also a gene given in transcript rates tables, but we'll drop it here
    melted_mutrates = melt.data.table(mutrates[relevant_regions, -"gene", on = "pid"], id.vars = c("pid"), variable.factor = F)
    setnames(melted_mutrates, "pid", "region")
  } else {
    melted_mutrates = melt.data.table(mutrates[relevant_regions, on = "gene"], id.vars = c("gene"), variable.factor = F)
    setnames(melted_mutrates, "gene", "region")
  }
  setnames(melted_mutrates, c("variable", "value"), c("gene_rate_grp", "raw_rate"))
  sample_region_rates[melted_mutrates, raw_rate := raw_rate, on = c('gene_rate_grp', 'region')]

  cds_trinuc_comp = get_ref_data(cesa, "gene_trinuc_comp")

  # Hash trinuc rates for faster runtime with large data sets (where there could be millions of queries of trinuc_mat)
  # Also creating a melted version for individual rate lookups.
  trinuc_rates = new.env(parent = emptyenv())
  trinuc_mat = cesa@trinucleotide_mutation_weights$trinuc_proportion_matrix[samples$Unique_Patient_Identifier, , drop = F]
  for (row in rownames(trinuc_mat)) { 
    trinuc_rates[[row]] = unname(trinuc_mat[row, ])
  }
  trinuc_by_sample = melt(as.data.table(trinuc_mat, keep.rownames = 'Unique_Patient_Identifier'),
                          id.vars = 'Unique_Patient_Identifier', variable.name = 'trinuc_mut')
  
  # Dot product of trinuc comp and patient's expected relative trinuc rates yields the denominator
  # for site-specific mutation rate calculation; numerator is raw gene rate multiplied by the patient's relative rate for the site's trinuc context
  # (this last value gets multipled in by get_baseline functions below)
  sample_region_rates[, aggregate_rate := raw_rate / sum(cds_trinuc_comp[[region]] * trinuc_rates[[Unique_Patient_Identifier]]), by = c("region", "Unique_Patient_Identifier")]
  setkey(sample_region_rates, "region")
  
  if(length(aac_ids > 0)) {
    mutations$aac_snv_key[mutations$snv, trinuc_mut := trinuc_mut, on = 'snv_id']
    trinuc_mut_by_aac = mutations$aac_snv_key[aac_ids, .(aac_id, trinuc_mut), on = 'aac_id']

    if (using_pid) {
      mutations$amino_acid_change[, region := pid]
    } else {
      mutations$amino_acid_change[, region := gene]
    }
    # for each sample, multiply agg_rate by relative trinuc rate of the context for all AAC sites
    if(length(aac_ids) > 1000) {
      message("Calculating baseline rates for amino-acid-changing mutations...")
    }

    trinuc_mut_by_sample_per_aac = merge.data.table(trinuc_mut_by_aac, trinuc_by_sample, by = 'trinuc_mut', allow.cartesian = T)
    sample_rates = trinuc_mut_by_sample_per_aac[, .(rate = sum(value)), by = c('aac_id', 'Unique_Patient_Identifier')]
    sample_rates[mutations$amino_acid_change, region := region, on = 'aac_id']
    sample_rates[sample_region_rates, final_rate := rate * aggregate_rate, on = c('region', 'Unique_Patient_Identifier')]
    final_aac_rates = dcast.data.table(sample_rates,  Unique_Patient_Identifier ~ aac_id, value.var = 'final_rate')
  } else {
    final_aac_rates = data.table(Unique_Patient_Identifier = samples$Unique_Patient_Identifier)
  }

  # repeat with SNVs (slightly different handling since SNVs can have more than one gene/pid associated)
  if (length(snv_ids) > 0) {
    if(length(snv_ids) > 1000) {
      message("Calculating baseline rates for noncoding SNVs...")
    }
    trinuc_mut_by_snv = mutations$snv[snv_ids, .(snv_id, trinuc_mut), on = 'snv_id']
    if (using_pid) {
      region_by_snv = mutations$snv[snv_ids, .(region = unlist(nearest_pid)), by = "snv_id"]
    } else {
      region_by_snv = mutations$snv[snv_ids, .(region = unlist(genes)), by = "snv_id"]
    }
    sample_rates_snv = merge.data.table(trinuc_mut_by_snv, trinuc_by_sample, by = 'trinuc_mut', allow.cartesian = T)
    
    
    rates_by_snv = merge.data.table(region_by_snv, sample_region_rates[, .(region, Unique_Patient_Identifier, aggregate_rate)],
                     by = 'region', all.x = T, all.y = F, allow.cartesian = T)
    
    # SNVs that cover multiple regions will get a rate averaged over those regions. In the future this may change.
    rates_by_snv = rates_by_snv[, .(snv_id, Unique_Patient_Identifier, aggregate_rate = mean(aggregate_rate)), by = c('snv_id', 'Unique_Patient_Identifier')]
    
    sample_rates_snv[rates_by_snv, final_rate := value * aggregate_rate, on = c('snv_id', 'Unique_Patient_Identifier')]
    final_snv_rates = dcast.data.table(sample_rates_snv,  Unique_Patient_Identifier ~ snv_id, value.var = 'final_rate')
  } else {
    final_snv_rates = data.table(Unique_Patient_Identifier = samples$Unique_Patient_Identifier)
  }
  
  # Combine all rates
  baseline_rates = merge.data.table(final_aac_rates, final_snv_rates, by = 'Unique_Patient_Identifier')
  setcolorder(baseline_rates, "Unique_Patient_Identifier")
  return(baseline_rates)
}


