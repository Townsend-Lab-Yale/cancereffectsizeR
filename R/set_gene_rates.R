#' Assign pre-calculated regional mutation rates
#'
#' This function allows you to specify regional rates of mutation--calculated
#' however you like--to samples in your CESAnalysis. Rates can be assigned to all samples
#' or to specified samples.
#' 
#' Provide rates in a data.table with two columns: gene name or protein ID (character) and
#' rate (numeric, non-negative). Gene names or protein IDs must match those in CESAnalysis
#' reference data. (Some reference data sets, such as ces.refset.hg19, only allow
#' gene-level rates.) If you don't want to supply rates for every gene, set
#' \code{missing_genes_take_nearest = T} to have each missing gene or coding region take
#' the rate of the nearest non-missing one.
#' 
#' 
#' @param cesa CESAnalysis object
#' @param rates A two-column data.table with either gene name or protein_id in column 1
#'   and rate in column 2
#' @param samples Which samples the input rates apply to. Defaults to all samples. Can be
#'   a vector of Unique_Patient_Identifiers, or a data.table containing rows from the
#'   CESAnalysis sample table.
#' @param missing_genes_take_nearest Set to TRUE to have each gene/protein_id missing from
#'   rates take the rate of the nearest non-missing gene/protein.
#' @export
set_gene_rates = function(cesa = NULL, rates = NULL, samples = character(), 
                          missing_genes_take_nearest = FALSE) {
  stopifnot(is(cesa, "CESAnalysis"),
            is(rates, "data.table"),
            is(missing_genes_take_nearest, "logical"),
            length(missing_genes_take_nearest) == 1)
  
  cesa = copy_cesa(cesa)
  if(cesa@samples[, .N] == 0) {
    stop("No MAF data loaded. Gene rates should be assigned after loading all variant data.")
  } 
  
  curr_sample_info = select_samples(cesa, samples)

  
  # Don't allow rates to be overwritten (too confusing for now to risk overlapping sample group specification in sequential calls)
  if (! all(is.na(curr_sample_info$gene_rate_grp))) {
    msg = paste0("Gene mutation rates have already been calculated/assigned for some or all samples specified. ",
         "Use clear_gene_rates() first if you want to assign new rates.")
    stop(pretty_message(msg, emit = F))
  }
  
  # validate rates
  rates = copy(rates)
  if(ncol(rates) != 2) {
    stop("rates should be a two-column table giving genes (or protein IDs) and their associated rates.")
  }
  if(! 'rate' %in% names(rates)) {
    stop('Expected a column named \"rate\" in input gene rate table.')
  }
  setnames(rates, setdiff(names(rates), 'rate'), 'region') # gene or pid becomes region
  if (! is(rates$region, "character") | ! is(rates$rate, "numeric")) {
    stop("In rates, gene (or protein_id) and rate should be character and numeric, respectively")
  }
  if(anyNA(rates)) {
    stop("There are NA entries in rates.")
  }
  if(any(duplicated(rates$region))) {
    stop("There are duplicate gene entries in rates.")
  }
  
  if(rates[, any(rate < 0)]) {
    stop("All rates must be non-negative.")
  }
  
  # Load refset genes and compare to rates
  gr_genes = get_ref_data(cesa, "gr_genes")
  
  
  # Presence of gene column (rather than just names) means that names gives pid
  pid_available = FALSE
  all_genes = unique(gr_genes$names)
  if ("gene" %in% names(GenomicRanges::mcols(gr_genes))) {
    pid_available = TRUE
    all_pid = all_genes
    all_genes = unique(gr_genes$gene)
  }
  
  using_pid = FALSE
  all_regions = all_genes
  if(length(setdiff(rates$region, all_genes)) > 0) {
    if (pid_available && length(setdiff(rates$region, all_pid)) == 0) {
        using_pid = TRUE
        all_regions = all_pid
    } else {
      stop("There are gene names in rates that are absent from the CESAnalysis's reference data set")
    }
  }
  
  missing_regions = setdiff(all_regions, rates$region)
  if (length(missing_regions) > 0) {
    if (missing_genes_take_nearest == FALSE) {
      stop("Not all reference genes (or protein_id) present in rates. To have each gene/pid take the rate of the nearest ",
           "non-missing one, re-run with missing_genes_take_nearest = TRUE")
    } else {
      tmp = as.data.table(gr_genes)
      # names gives the regions we're using, unless gr_genes contains protein_ids and we're only using genes
      if(pid_available) {
        if (using_pid) {
          tmp[, gene := NULL]
        } else {
          tmp[, names := NULL]
          setnames(tmp, 'gene', 'names')
        }
      }
      cds_intervals = tmp[, .(chr = as.character(seqnames[1]), start = min(start), end = max(end)), by = 'names']
      setnames(cds_intervals, "names", "region")
      cds_intervals[, center := start + trunc((end - start)/2)]
      cds_intervals = cds_intervals[order(chr, center)]
      present_regions = cds_intervals[! region %in% missing_regions]
      not_present = cds_intervals[region %in% missing_regions]
      
      # find nearest genes by chr/center to the missing genes from the present genes,
      # and eliminate the rare tie by taking the first match for each gene
      nearest = present_regions[not_present, on = c("chr", "center"), roll = "nearest"]
      nearest = nearest[! duplicated(i.region)]
      
      setkey(rates, "region")
      nearest_rates = rates[nearest$region, -"region"]
      missing_rates = cbind(region = nearest$i.region, nearest_rates)
      rates = rbind(rates, missing_rates)
    }
  }
  
  curr_rate_group = as.integer(max(c(0, na.omit(cesa@samples$gene_rate_grp))) + 1)
  rate_grp_colname = paste0('rate_grp_', curr_rate_group)
  
  if (using_pid) {
    setnames(rates, c("region", "rate"), c("pid", rate_grp_colname))
    pid_to_gene = unique(as.data.table(GenomicRanges::mcols(gr_genes)))
    setnames(pid_to_gene, 'names', 'pid')
    rates[pid_to_gene, gene := gene, on = 'pid']
    setcolorder(rates, c('pid', 'gene'))
  } else {
    setnames(rates, c("region", "rate"), c("gene", rate_grp_colname))
  }
  
  cesa@samples[curr_sample_info$Unique_Patient_Identifier, gene_rate_grp := curr_rate_group, on = "Unique_Patient_Identifier"]
  
  if (curr_rate_group > 1) {
    if (using_pid) {
      cesa@mutrates = cesa@mutrates[rates[, - "gene"], on = "pid"]
    } else {
      cesa@mutrates = cesa@mutrates[rates, on = "gene"]
    }
  } else {
    cesa@mutrates = rates
  }
  rate_cols = names(cesa@mutrates)[names(cesa@mutrates) %like% 'rate']
  setcolorder(cesa@mutrates, setdiff(names(cesa@mutrates), rate_cols))
  
  cesa = update_cesa_history(cesa, match.call())
  return(cesa)
}