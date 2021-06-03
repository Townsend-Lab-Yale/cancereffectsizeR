#' Load MAF somatic mutation data
#' 
#' Load MAF data from a text file or data table into your CESAnalysis. If column names
#' don't match MAF format specifications (Chromosome, Start_Position, etc., with
#' Tumor_Sample_Barcode used as the sample ID column), you can supply your own column
#' names. When your CESAnalysis has defined sample groups (see \code{?CESAnalysis}),
#' specify "group_col". By default, data is assumed to be derived from whole-exome
#' sequencing. Whole-genome data and targeted sequencing data are also supported when the
#' "coverage" option is specified. If the data you are loading is from a different genome
#' build than your CESAnalysis, you can use the "chain_file" option to supply a UCSC-style
#' chain file, and your MAF coordinates will be automatically converted with
#' rtracklayer's version of liftOver.
#' 
#' @param cesa the CESAnalysis object to load the data into
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.table or data.frame format
#' @param sample_col column name with sample ID data (Tumor_Sample_Barcode or Unique_Patient_Identifier)
#' @param chr_col column name with chromosome data  (Chromosome)           
#' @param start_col column name with start position (Start_Position)
#' @param ref_col column name with reference allele data (Reference_Allele)
#' @param tumor_allele_col column name with alternate allele data; by default,
#'   values from Tumor_Seq_Allele2 and Tumor_Seq_Allele1 columns are used.
#' @param group_col column in MAF with sample group labels (see \code{?CESAnalysis})
#' @param coverage exome, genome, or targeted (default exome)
#' @param covered_regions optional for exome, required for targeted: a GRanges object or a
#'   BED file of covered intervals matching the CESAnalysis genome
#' @param covered_regions_name a name describing the covered regions (e.g.,
#'   "my_custom_targeted_regions"); required when covered_regions are supplied
#' @param covered_regions_padding How many bases (default 0) to expand start and end of
#'   each covered_regions interval, to include variants called just outside of targeted
#'   regions. Consider setting from 0-100bp, or up to the sequencing read length. If the
#'   input data has been trimmed to the targeted regions, leave set to 0.
#' @param chain_file a LiftOver chain file (text format, name ends in .chain) to convert MAF
#'   records to the genome build used in the CESAnalysis.
#' @param enforce_default_exome_coverage When loading default exome data, exclude records
#'   that aren't covered in the default exome capture intervals included with
#'   CES genome reference data (default FALSE).
#' @return CESAnalysis with the specified MAF data loaded
#' @export
load_maf = function(cesa = NULL, maf = NULL, sample_col = "Tumor_Sample_Barcode", chr_col = "Chromosome", start_col = "Start_Position",
                    ref_col = "Reference_Allele", tumor_allele_col = "guess", coverage = "exome", covered_regions = NULL,
                    covered_regions_name = NULL, covered_regions_padding = 0, group_col = NULL, chain_file = NULL, enforce_default_exome_coverage = FALSE) {
  
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  
  if(cesa@advanced$locked) {
    stop("You can't load more MAF data since you've already calculated some mutation rates. Create a new one if necessary.", call. = F)
  }
  cesa = update_cesa_history(cesa, match.call())
  bsg = get_cesa_bsg(cesa)
  
  # validate chain_file (presence means liftOver must run)
  need_to_liftOver = FALSE
  if(! is.null(chain_file)) {
    need_to_liftOver = TRUE
    if(!is(chain_file, "character") || length(chain_file) != 1 || !(endsWith(chain_file, ".chain"))) {
      stop("Argument chain_file expected to be the filename/path of a text file ending in .chain")
    }
    if(! file.exists(chain_file)) {
      stop("The chain_file specified could not be found; check the file path.")
    }
  }
  
  if (is.null(maf)) {
    stop("Supply MAF data via maf=[file path or data.table/data.frame].")
  }
  
  if (is.null(group_col) & length(cesa@groups) != 1) {
    stop("You must specify group_col in the MAF since this CESAnalysis specifies sample groups.")
  }
  
  if (! is.null(group_col) & length(cesa@groups) == 1) {
    stop(paste0("This CESAnalysis does not specify sample groups, so you can't use the \"group_col\" argument.\n",
                "Create a new CESAnalysis with \"sample_groups\" specified to include this information."))
  }
  
  # give a warning if interval padding is really high
  if (covered_regions_padding > 1000) {
    warning(sprintf("covered_regions_padding=%d is awfully high!", covered_regions_padding), call. = F)
  }
  
  # Validate covered_regions
  previous_covered_regions_names = c(names(cesa@coverage$exome), names(cesa@coverage$targeted)) # may be NULL
  if (! is.character(coverage) || ! coverage %in% c("exome", "genome", "targeted") || length(coverage) > 1) {
    stop("Argument coverage must be \"exome\", \"genome\", or \"targeted\"")
  }
  
  if (! coverage %in% names(cesa@coverage)) {
    cesa@coverage[[coverage]] = list()
  }
  
  
  # validate covered_regions_name (supplied if and only if covered_regions is)
  if (! is.null(covered_regions) & is.null(covered_regions_name)) {
    stop("You must supply a name for your covered_regions using covered_regions_name = ...")
  }
  if (is.null(covered_regions) && ! is.null(covered_regions_name)) {
    stop("covered_regions_name was supplied, but covered_regions wasn't.")
  }
  if (! is.null(covered_regions_name)) {
    if (! is(covered_regions_name, "character") || length(covered_regions_name) > 1) {
      stop("covered_regions_name should be a 1-length character vector, just a name to use for your covered_regions")
    }
  }
  
  if(is.null(covered_regions) && (length(covered_regions_padding) != 1 || covered_regions_padding != 0)) {
    stop("You can't use covered_regions_padding without supplying covered_regions.", call. = F)
  }
  
  # validate covered_regions (required for targeted, optional for exome/genome)
  if (coverage == "genome" & is.null(covered_regions)) {
    covered_regions_name = "genome"
  }
  
  if (coverage == "targeted" & is.null(covered_regions) ) {
    stop("can't load targeted data without covered_regions (see docs)")
  }
  
  if (! is(enforce_default_exome_coverage, "logical") || length(enforce_default_exome_coverage) != 1) {
    stop("enforce_default_exome_coverage should be T/F")
  }
  
  if (coverage == "exome" & is.null(covered_regions)) {
    covered_regions_name = "exome"
    if(! check_for_ref_data(cesa, "default_exome_gr")) {
      stop("This genome has no default exome intervals, so to load exome data you must supply covered_regions (see docs)")
    }
  }
  
  # Determine whether we're using "exome" or "exome+" for default exome data
  # Usually, the lenient exome+ option is used, but the choice must be consistent throughout a CESAnalysis
  if(covered_regions_name == "exome") {
    # Whether or not CESAnalysis uses "exome+", once the first default exome data has been loaded,
    # there will always be default exome intervals under "exome" in the coverage list
    if (! "exome" %in% names(cesa@coverage[["exome"]])) {
      cesa@advanced$using_exome_plus = ! enforce_default_exome_coverage
      covered_regions = get_ref_data(cesa, "default_exome_gr")
      cesa@coverage$exome[["exome"]] = covered_regions
      
    } else if (enforce_default_exome_coverage == TRUE & cesa@advanced$using_exome_plus) {
      stop("You can't have enforce_default_exome_coverage = TRUE because you previously loaded default exome data\n",
           "with enforce_default_exome_coverage = FALSE")
    } else if (enforce_default_exome_coverage == FALSE & ! cesa@advanced$using_exome_plus) {
      stop("You can't have enforce_default_exome_coverage == FALSE because you previously loaded default exome data\n",
           "with enforce_default_exome_coverage = TRUE.")
    } else {
      covered_regions = cesa@coverage$exome[["exome"]]
    }
    
    # For exome+, use GRanges if available
    # Otherwise, for exome or exome+, load the default exome intervals
    if (! enforce_default_exome_coverage) {
      covered_regions_name = "exome+"
      if("exome+" %in% names(cesa@coverage$exome)) {
        covered_regions = cesa@coverage$exome[["exome+"]]
      } else {
        cesa@coverage[["exome"]][["exome+"]] = covered_regions # that is, the default exome regions
      }
    } else {
      if("exome" %in% names(cesa@coverage$exome)) {
        covered_regions = cesa@coverage$exome[["exome"]]
      } else {
        covered_regions = get_ref_data(cesa, "default_exome_gr")
        cesa@coverage$exome[["exome"]] = covered_regions
      }
    }
    pretty_message("Assuming this data has default exome coverage (it's better to supply covered intervals if you have them; see docs)...")
  } else if (! is.null(covered_regions)) {
    cesa = .add_covered_regions(cesa = cesa, covered_regions = covered_regions, covered_regions_padding = covered_regions_padding, 
                               coverage_type = coverage, covered_regions_name = covered_regions_name)
  }
  
  refset = .ces_ref_data[[cesa@ref_key]]
  read_args = list(maf = maf, refset_env = refset,
                   sample_col = sample_col, chr_col = chr_col, start_col = start_col,
                   ref_col = ref_col, tumor_allele_col = tumor_allele_col,
                   chain_file = chain_file)
  if (! is.null(group_col)) {
    read_args = c(read_args, list(more_cols = group_col))
  }
  maf = do.call(read_in_maf, args = read_args)
  

  # Set aside records with problems and notify user
  initial_num_records = maf[, .N]
  excluded = maf[! is.na(problem), .(Unique_Patient_Identifier, Chromosome, Start_Position, 
                                     Reference_Allele, Tumor_Allele, problem)]
  maf = maf[is.na(problem), -"problem"]
  maf = identify_maf_variants(maf) # add variant_type/variant_id columns

  num_excluded = excluded[, .N]
  if(num_excluded > 0) {
    msg = paste0(num_excluded, " of ", initial_num_records, " MAF records (", 
                 sprintf("%.1f", 100 * num_excluded / initial_num_records), '%) ',
                 "had problems and were excluded: ")
    problem_summary = excluded[, .(num_records = .N), by = "problem"]
    message(crayon::black(paste0(utils::capture.output(print(problem_summary, row.names = F)), collapse = "\n")))
    
    if(num_excluded / initial_num_records > .05) {
      warning("More than 5% of input records had problems.")
    }
    nt = c("A", "C", "G", "T")
    snv_mismatch = excluded[problem == "reference_mismatch" & Reference_Allele %in% nt & Tumor_Allele %in% nt, .N]
    if(snv_mismatch > 0) {
      msg = paste0(snv_mismatch, " SNV variants were excluded for having reference alleles that do not match the reference genome. You should probably figure out why",
                          " and make sure that the rest of your data set is okay to use before continuing.")
      warning(pretty_message(msg, emit = F))
    }
  }
  setnames(excluded, "problem", "Exclusion_Reason")
  
  # collect sample group information
  if (! is.null(group_col)) {
    if (is.factor(maf[[group_col]])) {
      warning("You supplied a sample group column as a factor, but it was converted to character.")
    }
    sample_groups = as.character(maf[[group_col]])
    if(any(is.na(sample_groups))) {
      stop("Error: There are NA values in your sample groups column.")
    }
    maf[, (group_col) := NULL]
  } else {
    sample_groups = cesa@groups[1] # indicates a stageless analysis
  }
  
  new_samples = data.table(Unique_Patient_Identifier = maf$Unique_Patient_Identifier, group = sample_groups)
  new_samples = new_samples[, .(group = unique(group)), by = "Unique_Patient_Identifier"]

  new_samples[, coverage := coverage]
  new_samples[, covered_regions := covered_regions_name]
  
  # ensure no sample has an illegal group
  bad_groups = setdiff(new_samples[, unique(group)], cesa@groups)
  if(length(bad_groups) > 0) {
    stop(paste0("The following groups were not declared in your CESAnalysis, but they were found in your MAF groups column:\n",
                paste(bad_groups, collapse = ", ")))
  }
  # see if any sample appears more than once in sample table (happens when one sample has multiple listed groups)
  repeated_samples = new_samples[duplicated(Unique_Patient_Identifier), unique(Unique_Patient_Identifier)]
  if(length(repeated_samples) > 0) {
    stop(paste0("The following samples are associated with multiple groups in the input data:\n", paste(repeated_samples, collapse=", ")))
  }
  
  # make sure no new samples were already in the previously loaded MAF data
  if(cesa@samples[, .N] > 0) {
    repeat_samples = intersect(cesa@samples[, Unique_Patient_Identifier], new_samples[, Unique_Patient_Identifier])
    if (length(repeat_samples) > 0) {
      stop(paste0("Error: Can't load MAF data because some sample IDs already appear in previously loaded data.\n",
                  "Either merge these data sets manually or remove duplicated samples: ",
                  paste(repeat_samples, collapse = ", ")))
    }
  }
  
  # notify the user if some of the declared groups don't appear in the data at all
  if(length(cesa@groups) > 1) {
    missing_groups = cesa@groups[! cesa@groups %in% new_samples[,unique(group)]]
    if (length(missing_groups) > 0) {
      msg = paste0("The following groups were declared in your CESAnalysis, but they weren't present in the MAF data: \n",
                     paste(missing_groups, collapse = ", "))
      pretty_message(msg)
    }    
  }

  # remove any MAF records that are not in the coverage, unless default exome with enforce_default_exome_coverage = FALSE
  if (covered_regions_name == "genome") {
    num_uncovered = 0
  } else {
    maf_grange = GenomicRanges::makeGRangesFromDataFrame(maf, seqnames.field = "Chromosome", start.field = "Start_Position", 
                                                         end.field = "Start_Position")
    
    # In an "exome+" data set, compare coverage to default exome instead of whatever the current exome+ intervals are
    if (covered_regions_name == "exome+") {
      # equivalent to %within%, but avoids importing
      is_uncovered = ! IRanges::overlapsAny(maf_grange, cesa@coverage[["exome"]][["exome"]], type = "within")
    } else {
      is_uncovered = ! IRanges::overlapsAny(maf_grange, cesa@coverage[[coverage]][[covered_regions_name]], type = "within")
    }
    num_uncovered = sum(is_uncovered)
  }

  if (num_uncovered > 0) {
    total = nrow(maf)
    percent = round((num_uncovered / total) * 100, 1)
    if (covered_regions_name == "exome+") {
      # merge previous exome+ covered_regions with the coverage of the new data
      covered_regions =  GenomicRanges::reduce(GenomicRanges::union(covered_regions, maf_grange[is_uncovered]))
      cesa@coverage[[coverage]][["exome+"]] = covered_regions

      # warn if a lot of records are in uncovered areas; may indicate whole-genome data or low quality exome data
      if (percent > 10) {
        warning(paste0("More than 10% of MAF records are not within the CES genome's default exome intervals.\n",
                       "Could this be whole-genome data? Or if you know the true covered regions, supply them\n",
                       "with the covered_regions argument."))
      }
      msg = paste0("Note: ", format(num_uncovered, big.mark = ','), " MAF records (", percent, 
                   "%) are outside the CESAnalysis's default exome definitions; expanded exome intervals to include them.")
      pretty_message(msg)
    } else {
      uncovered = maf[is_uncovered, -c("variant_type", "variant_id")]
      uncovered$Exclusion_Reason = paste0("uncovered_in_", covered_regions_name)
      maf = maf[!is_uncovered]
      excluded = rbind(excluded, uncovered)
      message(paste0("Note: ", num_uncovered, " MAF records out of ", total, " (", percent, 
                     "%) are at loci not covered in the input covered_regions, so they have been excluded."))
    }
  }
  
  # drop any samples that had all mutations excluded
  new_samples = new_samples[Unique_Patient_Identifier %in% maf$Unique_Patient_Identifier]
  cesa@samples = rbind(cesa@samples, new_samples)
  setcolorder(cesa@samples, c("Unique_Patient_Identifier", "coverage", "covered_regions", "group"))
  setkey(cesa@samples, "Unique_Patient_Identifier")
  
  if (nrow(excluded) > 0) {
    colnames(excluded) = c(colnames(maf)[1:5], "Exclusion_Reason")
    cesa@excluded = rbind(cesa@excluded, excluded) 
  }
  
  # Set aside new variants for annotation (notably, before MNV prediction; we'll still annotate those as SNVs)
  # To-do: also leave out indels that have previously been annotated (okay to re-annotate for now)
  to_annotate = maf[! cesa@mutations$snv$snv_id, on = 'variant_id']
  message("Annotating variants...")
  annotations = annotate_variants(refset = refset, variants = to_annotate)
  
  # Check annotations for any SNVs with ambiguous trinuc context; these must be set aside
  snv_table = annotations$snv
  aac_table = annotations$amino_acid_change
  bad_trinuc_context = which(is.na(snv_table$trinuc_mut))
  num_bad = length(bad_trinuc_context)
  if (num_bad > 0) {
    bad_ids = snv_table[bad_trinuc_context, snv_id]
    bad_trinuc_context_maf = maf[bad_ids, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele), on = 'variant_id']
    maf = maf[! bad_ids, on = 'variant_id']
    msg = paste0("Note: ", num_bad, " MAF records excluded due to ambiguous trinucleotide context ",
                  "(likely N's in the reference genome).")
    pretty_message(msg)
    bad_trinuc_context_maf$Exclusion_Reason = "ambiguous_trinuc_context"
    cesa@excluded = rbind(cesa@excluded, bad_trinuc_context_maf)
    
    # for simplicity, remove the bad record from SNV and AAC tables
    bad_aa = unlist(snv_table[bad_trinuc_context, assoc_aac])
    snv_table = snv_table[! bad_trinuc_context]
    if(aac_table[, .N] > 0 & length(bad_aa) > 0) {
      aac_table = aac_table[! bad_aa]
    }
  }
  cesa@mutations[["amino_acid_change"]] = unique(rbind(cesa@mutations$amino_acid_change, aac_table, fill = T), by = "aac_id")
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  cesa@mutations[["snv"]] = unique(rbind(cesa@mutations$snv, snv_table, fill = T), by = "snv_id")
  setkey(cesa@mutations$snv, "snv_id")
  
  # add genes and assoc_aac to MAF records (may stop including these in near future)
  # use of _tmp names required as of data.table 1.13.2 to keep join from failing
  column_order = copy(names(maf))
  maf[, c("genes_tmp", "assoc_aac_tmp") := list(list(NA_character_), list(NA_character_))]
  maf[cesa@mutations$snv, c("genes_tmp", "assoc_aac_tmp") := list(genes, assoc_aac), on = c(variant_id = "snv_id")]
  
  # No gene for intergenic SNVs (SNV annotation table gives nearest gene, but we won't show that in MAF table)
  intergenic_snv = cesa@mutations$snv[intergenic == T, snv_id]
  maf[intergenic_snv, genes_tmp := list(NA_character_), on = "variant_id"]
  
  
  # temporary way of annotating non-SNVs
  non_snv = maf[variant_type != 'snv']
  if (non_snv[, .N] > 0) {
    grt = as.data.table(refset$gr_genes) # Name will change to gr_cds later
    if ("gene" %in% names(grt)) {
      grt[, names := gene] # use gene field instead of names field (applies when CDS gr has multiple CDS per gene)
      grt[, gene := NULL]
    }
    setnames(grt, c("seqnames", "start", "end", "names"), 
             c("Chromosome", "Start_Position", "End_Position", "gene"))
    setkey(grt, Chromosome, Start_Position, End_Position)
    non_snv[, End_Position := Start_Position] # okay for now
    non_snv_genes = foverlaps(non_snv, grt)
    # use Start_Position from MAF, not gr
    non_snv_genes[, Start_Position := i.Start_Position]
    non_snv_genes = non_snv_genes[, .(genes = list(unique(gene))), by = c("Chromosome", "Start_Position")]
    non_snv_genes[, is_snv := F]
    maf[, is_snv := variant_type == 'snv']
    maf[non_snv_genes, c("genes_tmp", "assoc_aac_tmp") := .(genes, list(NA_character_)),
        on = c("is_snv", "Chromosome", "Start_Position")]
    maf[, is_snv := NULL]
  }
  
  setnames(maf, c("genes_tmp", "assoc_aac_tmp"), c("genes", "assoc_aac"))
  setcolorder(maf, column_order) # put original columns back in the front
  
  # Same-sample variants with 2bp of other variants get set aside as likely MNVs
  # MNVs are only possible in sample/chromosome combinations with more than one MAF record
  poss_mnv = maf[order(Unique_Patient_Identifier, Chromosome, Start_Position)][, .SD[.N > 1], 
                                                                                    by = c("Unique_Patient_Identifier", "Chromosome")]
  if (poss_mnv[, .N] > 0) {
    poss_mnv[, dist_to_prev := c(Inf, diff(Start_Position)), by = c("Unique_Patient_Identifier", "Chromosome")]
    poss_mnv[, dist_to_next := c(dist_to_prev[2:.N], Inf), by = c("Unique_Patient_Identifier", "Chromosome")]
    poss_mnv[dist_to_prev < 3 | dist_to_next < 3, is_mnv := T]
    
    # organize into groups of likely multi-nucleotide events
    mnv = poss_mnv[is_mnv == T]
    mnv[, start_of_group := dist_to_prev > 2]
    mnv[, mnv_group := cumsum(start_of_group)]
    
    # Groups of 2 consecutive SNVs are double-base substitutions
    # Create DBS entries suitable for MAF table
    mnv[, is_dbs := .N == 2 && diff(Start_Position) == 1 && variant_type[1] == 'snv' && variant_type[2] == 'snv', by = mnv_group]
    dbs = mnv[is_dbs == T, 
                  .(Unique_Patient_Identifier = Unique_Patient_Identifier[1], 
                    Chromosome = Chromosome[1], Start_Position = Start_Position[1],
                    Reference_Allele = paste(Reference_Allele, collapse = ''),
                    Tumor_Allele = paste0(Tumor_Allele, collapse = ''), variant_type = 'dbs',
                    v1 = variant_id[1], v2 = variant_id[2]), by = mnv_group][, -"mnv_group"]
    dbs[, dbs_id := paste0(Chromosome, ':', Start_Position, '_', Reference_Allele, '>', Tumor_Allele)]
    dbs[maf, genes := genes, on = c(v1 = 'variant_id')]
    dbs[, assoc_aac := list(NA_character_)]
    
    # Remove the SNV entries that form the new DBS variants
    maf[dbs, dbs_id := i.dbs_id, on = c("Unique_Patient_Identifier", variant_id = "v1")]
    maf[dbs, dbs_id := i.dbs_id, on = c("Unique_Patient_Identifier", variant_id = "v2")]
    maf_dbs_ind = maf[! is.na(dbs_id), which = T]
    maf_dbs = maf[maf_dbs_ind, .(Unique_Patient_Identifier, Chromosome, Start_Position, 
                                 Reference_Allele, Tumor_Allele, 
                                 Exclusion_Reason = paste0("replaced_with_", dbs_id))]
    
    if (all(c("prelift_chr", "prelift_start", "liftover_strand_flip") %in% names(maf))) {
      dbs[maf, c("prelift_chr", "prelift_start", "liftover_strand_flip") := list(prelift_chr, prelift_start, liftover_strand_flip), on = c(v1 = 'variant_id')]
    }
    
    maf = maf[! maf_dbs_ind, -"dbs_id"]

    # Add new DBS entries
    dbs[, c("v1", "v2") := NULL]
    setnames(dbs, "dbs_id", "variant_id")
    maf = rbind(maf, dbs)
    
    # For non-DBS mnvs, reclassify as "other"
    mnv = mnv[is_dbs == F]
    maf[mnv, variant_type := 'other', on = c('Unique_Patient_Identifier', 'variant_id')]
    num_new_mnv = uniqueN(mnv$mnv_group) # the DBS variants have already been removed
    cesa@excluded = rbind(cesa@excluded, maf_dbs)
    
    num_dbs = dbs[, .N]
    if(num_dbs > 0) {
      grammar = ifelse(num_dbs == 1, 'has', 'have')
      msg = paste0('Note: ', num_dbs, ' adjacent pairs of SNVs ', grammar, ' been re-annotated as doublet base substitutions (dbs).')
      pretty_message(msg)
    }
    if (num_new_mnv > 0) {
      msg = ifelse(num_dbs > 0, 'Additionally, ', 'Note: ')
      grammar = ifelse(num_new_mnv > 1, 's', '')
      msg = paste0(msg, num_new_mnv, ' group', grammar, ' of same-sample variants within 2 bp of each other have been reclassified as ',
              'combined variants of type \"other\" (because they probably did not occur independently).')
      pretty_message(msg)
    }
    maf[variant_type == 'other', c("variant_id", "assoc_aac") := list(NA_character_, list(NA_character_))] 
  }
  
  
  cesa@maf = rbind(cesa@maf, maf, fill = T)
  
  
  # Update coverage fields of annotation tables
  # Internal note: Confusingly, update_covered_in also has a side effect of updating
  # cached output of select_variants; this should change in the future.
  cesa = update_covered_in(cesa)

  current_snv_stats = maf[variant_type == "snv", .(num_samples = uniqueN(Unique_Patient_Identifier), num_snv = .N)]
  msg = paste0("Loaded ", format(current_snv_stats$num_snv, big.mark = ','), " SNVs from ", 
                 format(current_snv_stats$num_samples, big.mark = ','), " samples into CESAnalysis.")
  pretty_message(msg)
  
  return(cesa)
}



