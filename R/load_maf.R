#' Load MAF somatic mutation data
#' 
#' Load MAF data from a text file or data table into your CESAnalysis. Column names are
#' expected to match MAF format specifications (Chromosome, Start_Position, etc.). It's
#' recommended to use preload_maf() to prep the input (including, optionally, liftOver
#' conversion of genomic coordinates), but if you have clean MAF data, you can run this
#' function directly. By default, data is assumed to be derived from whole-exome
#' sequencing. Whole-genome data and targeted sequencing data are also supported when the
#' \code{coverage} option is specified.
#' 
#' @param cesa CESAnalysis.
#' @param maf Path of tab-delimited text file in MAF format, or an MAF in data.table or
#'   data.frame format.
#' @param maf_name Optionally, a name to identify samples coming from the current MAF. Used to
#'   populate the maf_source field of the CESAnalysis samples table.
#' @param sample_data_cols MAF columns containing sample-level data (e.g., tumor grade)
#'   that you would like to have copied into the CESAnalysis samples table.
#' @param coverage exome, genome, or targeted (default exome).
#' @param covered_regions optional for exome, required for targeted: a GRanges object or a
#'   BED file of covered intervals matching the CESAnalysis genome.
#' @param covered_regions_name a name describing the covered regions (e.g.,
#'   "my_custom_targeted_regions"); required when covered_regions are supplied.
#' @param covered_regions_padding How many bases (default 0) to expand start and end of
#'   each covered_regions interval, to include variants called just outside of targeted
#'   regions. Consider setting from 0-100bp, or up to the sequencing read length. If the
#'   input data has been trimmed to the targeted regions, leave set to 0.
#' @param enforce_default_exome_coverage When loading default exome data, exclude records
#'   that aren't covered in the default exome capture intervals included with CES genome
#'   reference data (default FALSE).
#' @return CESAnalysis with the specified MAF data loaded. The MAF data table includes
#'   CES-generated variant IDs, a list of all genes overlapping the site, and top_gene and
#'   top_consequence columns that give the most significant annotated coding changes for
#'   each mutation record. Annotation precedence is determined by MAF prevalence (usually
#'   equal), essential splice status, premature stop codon, nonsilent status, MAF mutation
#'   prevalence across the transcript (often favors longer transcripts), and finally
#'   alphabetical order. The columns are recalculated when more data is loaded, so changes
#'   in MAF prevalence can change which variants are highlighted. Note that
#'   \code{[CESAnalysis]$variants} contains more information about all top_consequence
#'   variants and all noncoding variants from the MAF.
#' @export
load_maf = function(cesa = NULL, maf = NULL, maf_name = character(), coverage = "exome", covered_regions = NULL,
                    covered_regions_name = NULL, covered_regions_padding = 0,
                    sample_data_cols = character(), enforce_default_exome_coverage = FALSE) {
  
  if (! is(cesa, "CESAnalysis")) {
    stop("cesa should be a CESAnalysis")
  }
  
  cesa = copy_cesa(cesa)
  cesa = update_cesa_history(cesa, match.call())
  bsg = get_cesa_bsg(cesa)
  
  if (is.null(maf)) {
    stop("Supply MAF data via maf=[file path or data.table/data.frame].")
  }
  
  if(! is.character(sample_data_cols)) {
    stop("sample_data_cols should be character if used (a vector of column names)")
  }
  
  # Validate sample_data_cols
  if (length(sample_data_cols) > 0) {
    sample_data_cols = unique(na.omit(sample_data_cols))
  } 
  
  
  # give a warning if interval padding is really high
  if (covered_regions_padding > 1000) {
    warning(sprintf("covered_regions_padding=%d is awfully high!", covered_regions_padding), call. = F)
  }
  
  # Validate covered_regions
  previous_covered_regions_names = c(names(cesa@coverage$exome), names(cesa@coverage$targeted)) # may be NULL
  if (! is.character(coverage) || length(coverage) > 1 || ! coverage %in% c("exome", "genome", "targeted")) {
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
      
      # Tolerate reference assembly increment (e.g., GRCh38p.13->GRCh37p.14)
      GenomeInfoDb::genome(covered_regions) = GenomeInfoDb::genome(bsg)[1] 
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
    pretty_message("Assuming this data has default exome coverage....")
  } else if (! is.null(covered_regions)) {
    cesa = .add_covered_regions(cesa = cesa, covered_regions = covered_regions, covered_regions_padding = covered_regions_padding, 
                               coverage_type = coverage, covered_regions_name = covered_regions_name)
  }
  
  refset = .ces_ref_data[[cesa@ref_key]]
  
  read_args = list(maf = maf, refset_env = refset,
                   sample_col = 'Unique_Patient_Identifier', chr_col = 'Chromosome', start_col = 'Start_Position',
                   ref_col = 'Reference_Allele', tumor_allele_col = 'guess', separate_old_problems = TRUE)

  if(length(sample_data_cols) > 0) {
    read_args = c(read_args, list(more_cols = sample_data_cols))
  }
  
  reserved_cols = c("Unique_Patient_Identifier", "Chromosome", "Start_Position", "Reference_Allele", 
                    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Allele", "Tumor_Sample_Barcode")
  illegal_sample_cols = intersect(reserved_cols, read_args[["more_cols"]])
  if(length(illegal_sample_cols) > 0) {
    stop("Can't use these column(s) as sample-level data columns: ", paste(illegal_sample_cols, collapse = ", "), '.')
  }

  if(! is.character(maf_name) || length(maf_name) > 1) {
    stop("maf_name should be 1-length character.")
  }
  if(length(maf_name) == 0) {
    maf_name = as.character(uniqueN(cesa@samples$maf_source) + 1)
  } else {
    if(! maf_name %ilike% '^[a-z0-9][-0-9a-z\\_\\.,\\+]*$') {
      stop('Invalid maf_name. You can use alphanumerics and ",", "+", "_", "-", ".".)')
    }
  }
  
  maf = do.call(read_in_maf, args = read_args)
  old_problems = data.table()
  if('old_problem' %in% names(maf)) {
    old_problems = maf[! is.na(old_problem), -"problem"]
    setnames(old_problems, 'old_problem', 'reason') # to put in exclusion table
    if(old_problems[, .N] > 0) {
      pretty_message("Found a \"problem\" column from preload_maf() and excluded records with problems.")
    }
    maf = maf[is.na(old_problem), -"old_problem"]
  }
  
  sample_data = NULL
  if(length(sample_data_cols) > 0) {
    missing_sample_data = setdiff(sample_data_cols, names(maf))
    if(length(missing_sample_data) > 0) {
      stop("Some sample_data_cols not present in MAF: ", paste(missing_sample_data, collapse = ', '), '.')
    }
    
    bad_cols = character()
    for (col in sample_data_cols) {
      col_good = maf[, .(good = uniqueN(.SD[[1]]) == 1), by = "Unique_Patient_Identifier", .SDcols = col][, all(good)]
      if (! col_good) {
        bad_cols = c(bad_cols, col)
      }
    }
    if(length(bad_cols) > 0) {
      stop("Sample data column(s) do not have consistent values within each sample: ", paste(bad_cols, collapse = ", "), '.')
    }
    
    # save sample data for later (only need 1 row per sample, since it's all sample-level data)
    sample_data = unique(maf[, .SD, .SDcols = c("Unique_Patient_Identifier", sample_data_cols)], by = "Unique_Patient_Identifier")
    maf = maf[, .SD, .SDcols = setdiff(names(maf), sample_data_cols)]
  }
  
  
  # Set aside records with problems and notify user
  initial_num_records = maf[, .N]
  excluded = maf[! is.na(problem), .(Unique_Patient_Identifier, Chromosome, Start_Position, 
                                     Reference_Allele, Tumor_Allele, variant_id, variant_type, problem)]
  maf = maf[is.na(problem), -"problem"]

  num_excluded = excluded[, .N]
  if(num_excluded > 0) {
    msg = paste0(num_excluded, " of ", initial_num_records, " MAF records (", 
                 sprintf("%.1f", 100 * num_excluded / initial_num_records), '%) ',
                 "had problems and were excluded: ")
    problem_summary = excluded[, .(num_records = .N), by = "problem"]
    if(Sys.getenv("RSTUDIO") == "1" && rstudioapi::isAvailable() &&rstudioapi::getThemeInfo()$dark) {
      message(crayon::white(paste0(utils::capture.output(print(problem_summary, row.names = F)), collapse = "\n")))
    } else {
      message(crayon::black(paste0(utils::capture.output(print(problem_summary, row.names = F)), collapse = "\n")))
    }
    if(num_excluded / initial_num_records > .05) {
      warning("More than 5% of input records had problems.")
    }
    nt = c("A", "C", "G", "T")
    sbs_mismatch = excluded[problem == "reference_mismatch" & Reference_Allele %in% nt & Tumor_Allele %in% nt, .N]
    if(sbs_mismatch / initial_num_records > .01) {
      msg = paste0(sbs_mismatch, " sbs variants were excluded for having reference alleles that do not match the reference genome. You should probably figure out why",
                          " and make sure that the rest of your data set is okay to use before continuing.")
      warning(pretty_message(msg, emit = F))
    }
  }
  setnames(excluded, "problem", "reason")
  
  if(old_problems[, .N] > 0) {
    excluded_cols = names(excluded)
    excluded = rbind(excluded, old_problems[, ..excluded_cols])
  }
  
  new_samples = data.table(Unique_Patient_Identifier = unique(maf$Unique_Patient_Identifier))
  new_samples[, coverage := coverage]
  new_samples[, covered_regions := covered_regions_name]
  
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
      if (percent > 20) {
        warning(paste0("More than 20% of MAF records are not within the CES genome's default exome intervals.\n",
                       "Could this be whole-genome data? Or if you know the true covered regions, supply them\n",
                       "with the covered_regions argument."))
      }
      msg = paste0("Note: ", format(num_uncovered, big.mark = ','), " MAF records (", percent, 
                   "%) are outside the CESAnalysis's default exome definitions; expanded exome intervals to include them.")
      pretty_message(msg)
    } else {
      uncovered = maf[is_uncovered]
      uncovered$reason = 'uncovered'
      maf = maf[!is_uncovered]
      excluded = rbind(excluded, uncovered[, names(excluded), with = F])
      pretty_message(paste0("Note: ", num_uncovered, " MAF records out of ", total, " (", percent, 
                     "%) are at loci not covered in the input covered_regions, so they have been excluded."))
    }
  }
  
  # merge in any sample-level data
  if (! is.null(sample_data)) {
    new_samples = merge.data.table(new_samples, sample_data, by = "Unique_Patient_Identifier")
  }
  new_samples[, maf_source := maf_name]
  
  cesa@samples = rbind(cesa@samples, new_samples, fill = TRUE)
  setcolorder(cesa@samples, c("Unique_Patient_Identifier", "coverage", "covered_regions"))
  setkey(cesa@samples, "Unique_Patient_Identifier")
  
  if (excluded[, .N] > 0) {
    cesa@excluded = rbind(cesa@excluded, excluded, fill = T) 
  }
  
  # Set aside new variants for annotation (notably, before MNV prediction; we'll still annotate those as sbs)
  # To-do: also leave out indels that have previously been annotated (okay to re-annotate for now)
  to_annotate = maf[! c(cesa@mutations$sbs$sbs_id, cesa@mutations$dbs$dbs_id), on = 'variant_id']
  pretty_message("Annotating variants...")
  annotations = annotate_variants(refset = refset, variants = to_annotate)
  
  # Check annotations for any sbs with ambiguous trinuc context; these must be set aside
  sbs_table = annotations$sbs
  aac_table = annotations$amino_acid_change
  bad_trinuc_context = which(is.na(sbs_table$trinuc_mut))
  num_bad = length(bad_trinuc_context)
  if (num_bad > 0) {
    bad_ids = sbs_table[bad_trinuc_context, sbs_id]
    bad_trinuc_context_maf = maf[bad_ids, .(Unique_Patient_Identifier, Chromosome, Start_Position, Reference_Allele, Tumor_Allele, variant_id, variant_type), on = 'variant_id']
    maf = maf[! bad_ids, on = 'variant_id']
    msg = paste0("Note: ", num_bad, " MAF records excluded due to ambiguous trinucleotide context ",
                  "(likely N's in the reference genome).")
    pretty_message(msg)
    bad_trinuc_context_maf$reason = "ambiguous_trinuc_context"
    cesa@excluded = rbind(cesa@excluded, bad_trinuc_context_maf)
    
    # For simplicity, remove the bad record from sbs and AAC tables
    bad_aa = annotations$aac_sbs_key[bad_ids, aac_id, on = 'sbs_id', nomatch = NULL]
    sbs_table = sbs_table[! bad_trinuc_context]
    if(aac_table[, .N] > 0 & length(bad_aa) > 0) {
      aac_table = aac_table[! bad_aa]
    }
  }
  cesa@mutations[["amino_acid_change"]] = unique(rbind(cesa@mutations$amino_acid_change, aac_table, fill = T), by = "aac_id")
  setkey(cesa@mutations$amino_acid_change, "aac_id")
  cesa@mutations[["sbs"]] = unique(rbind(cesa@mutations$sbs, sbs_table, fill = T), by = "sbs_id")
  setkey(cesa@mutations$sbs, "sbs_id")
  cesa@mutations[["aac_sbs_key"]] = unique(rbind(cesa@mutations$aac_sbs_key, annotations$aac_sbs_key))
  setkey(cesa@mutations$aac_sbs_key, 'aac_id')
  
  cesa@mutations[["dbs_codon_change"]] = unique(rbind(cesa@mutations$dbs_codon_change, annotations$dbs_codon_change), by = 'dbs_aac_id')
  cesa@mutations[["dbs"]] = unique(rbind(cesa@mutations$dbs, annotations$dbs), by = 'dbs_id')
  cesa@mutations[["aac_dbs_key"]] = unique(rbind(cesa@mutations$aac_dbs_key, annotations$aac_dbs_key))
  
  # add genes list to MAF records (may stop including in near future)
  # use of _tmp names required as of data.table 1.13.2 to keep join from failing
  column_order = copy(names(maf))
  maf[, genes_tmp := list(list(NA_character_))]
  maf[cesa@mutations$sbs,  genes_tmp := list(genes), on = c(variant_id = "sbs_id")]
  
  # No gene for intergenic sbs (sbs annotation table gives nearest gene, but we won't show that in MAF table)
  intergenic_sbs = cesa@mutations$sbs[intergenic == T, sbs_id]
  maf[intergenic_sbs, genes_tmp := list(NA_character_), on = "variant_id"]
  
  
  # temporary way of annotating non-sbs/DBS
  non_sbs_dbs = maf[variant_type != 'sbs']
  if (non_sbs_dbs[, .N] > 0) {
    grt = as.data.table(refset$gr_genes) # Name will change to gr_cds later
    if ("gene" %in% names(grt)) {
      grt[, names := gene] # use gene field instead of names field (applies when CDS gr has multiple CDS per gene)
      grt[, gene := NULL]
    }
    setnames(grt, c("seqnames", "start", "end", "names"), 
             c("Chromosome", "Start_Position", "End_Position", "gene"))
    setkey(grt, Chromosome, Start_Position, End_Position)
    non_sbs_dbs[, End_Position := Start_Position] # okay for now
    non_sbs_dbs_genes = foverlaps(non_sbs_dbs, grt)
    # use Start_Position from MAF, not gr
    non_sbs_dbs_genes[, Start_Position := i.Start_Position]
    non_sbs_dbs_genes = non_sbs_dbs_genes[, .(genes = list(unique(gene))), by = c("Chromosome", "Start_Position")]
    non_sbs_dbs_genes[, is_sbs := F]
    maf[, is_sbs := variant_type == 'sbs']
    maf[non_sbs_dbs_genes, genes_tmp := genes,
        on = c("is_sbs", "Chromosome", "Start_Position")]
    maf[, is_sbs := NULL]
  }
  
  setnames(maf, "genes_tmp", "genes")
  setcolorder(maf, column_order) # put original columns back in the front
  
  
  # Same-sample variants with 2bp of other variants get set aside as likely MNVs
  # MNVs are only possible in sample/chromosome combinations with more than one MAF record
  if(maf[, sum(variant_type == 'dbs')] == 0) {
    hidden_mnv = detect_mnv(maf)
    num_mnv = uniqueN(hidden_mnv$mnv_group)
    if(num_mnv > 0) {
      msg = paste0("There are ", num_mnv, " groups of same-sample variants within 2 bp of each other. ",
                   "These may not have resulted from independent events. Consider using preload_maf() ",
                   "to reclassify these variants as doublet substitutions or other multinucleotide variants.")
      warning(pretty_message(msg, emit = F))
    }
  }
  cesa@maf = rbind(cesa@maf, maf, fill = T)
  
  # Update coverage fields of annotation tables
  # Note that update_covered_in also updates cached variants table,
  # since the table includes some coverage info for convenience.
  cesa = update_covered_in(cesa)
  
  # Cached variants is a table of non-overlapping mutations (in terms of genomic position),
  # with only the "top" (tiebreaker-winning) variant at each site.
  # For sites that have coding effects, we will add this effect to the MAF table.
  # Edge case: No passing variants in MAF table means nothing to do. Also possible
  # to have no cached variants if all variants are of non-annotated type (e.g., "other").
  if(cesa@maf[, .N] > 0 && ! is.null(cesa@advanced$cached_variants)) {
    consequences = cesa@mutations$aac_sbs_key[cesa@advanced$cached_variants[variant_type == 'aac'], .(variant_name, sbs_id, gene), on = c(aac_id = 'variant_id')]
    cesa@maf[consequences, c("top_gene", "top_consequence") := list(gene, variant_name), on = c(variant_id = 'sbs_id')]
  }

  msg = paste0("Loaded ", format(maf[, .N], big.mark = ','), " variant records from ", 
                 format(new_samples[, .N], big.mark = ','), " samples into CESAnalysis.")
  pretty_message(msg)
  
  return(cesa)
}



