#' Make a PathScore input file from MAF data
#' 
#' Produce a file from MAF data for use with PathScore, a web tool for identifying significantly
#' altered pathways in cancer. See https://pathscore.publichealth.yale.edu/ and
#' https://doi.org/10.1093/bioinformatics/btw512.
#' 
#' @param maf An MAF-like data.table, as from preload_maf(). If the MAF has a column named \code{annot},
#'   this column will be preserved in output (since PathScore supports an optional annotation column
#'   with this name).
#' @param file Name/path where PathScore-formatted data should be written. Will be a tab-delimited text file.
#' @param genome The genome build associated with the MAF file. Must be "hg38" (default) or "hg19".
#' @return A copy of the PathScore-formatted data as a data.table.
#' @export
make_PathScore_input = function(maf, file = NULL, genome = 'hg38') {
  if(is.null(file)) {
    message("No filepath was specified (file argument NULL), so no PathScore input file will be created.")
  } else {
    if(! is(file, "character") || length(file) != 1) {
      stop("file should be a valid file path (1-length character) or left NULL.")
    }
    if(! grepl('\\.(txt|tsv)$', file)) {
      stop("filename should end in .txt or .tsv (it's going to be a tab-delimited text file).")
    }
    
    if (file.exists(file)) {
      stop(pretty_message(paste0('Specified output file ', file, ' already exists.'), emit = F))
    }
    dir = dirname(file)
    if(! dir.exists(dir)) {
      stop("The directory specified in the output file path does not exist.")
    }
    
    if(file.access(dir, 2) != 0) {
      if (dir == '.') {
        msg = paste0("You don't have write permissions in your working directory.",
                     " Change directories or specify a different path for your output filename.")
        stop(pretty_message(msg, emit = F))
      } else {
        stop("The directory specified in the output file path is not writeable.")
      }
    }
  }
  
  if(! is(maf, 'data.table')) {
    stop('maf should be type data.table')
  }
  maf = copy(maf)
  
  
  if (! 'Tumor_Allele' %in% names(maf) && 'Tumor_Seq_Allele2' %in% names(maf)) {
    maf[, Tumor_Allele := Tumor_Seq_Allele2]
  }
  
  if(! 'patient_id' %in% names(maf)) {
    if('Unique_Patient_Identifier' %in% names(maf)) {
      setnames(maf, 'Unique_Patient_Identifier', 'patient_id')
    } else if('Tumor_Sample_Barcode' %in% names(maf)) {
      setnames(maf, 'Tumor_Sample_Barcode', 'patient_id')
      message("Using column Tumor_Sample_Barcode as patient_id.")
    } 
  }
  
  required_maf_cols = c('Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Allele', 'patient_id')
  missing_cols = setdiff(required_maf_cols, names(maf))
  if(length(missing_cols) > 0) {
    missing_cols[missing_cols == 'patient_id'] = 'patient_id (or Unique_Patient_Identifier/Tumor_Sample_Barcode)'
    missing_cols[missing_cols == 'Tumor_Allele'] = 'Tumor_Allele (or Tumor_Seq_Allele2)'
    msg = paste0('Required MAF columns missing: ', paste(missing_cols, collapse = ', '), '.')
    stop(pretty_message(msg, emit = F))
  }
  
  maf = identify_maf_variants(maf)
  num_illegal = maf[variant_type == 'illegal', .N]
  if(num_illegal > 0) {
    warning(num_illegal, " MAF records appear invalid. These records will be excluded from output.")
  }
  maf = maf[variant_type != 'illegal']
  maf[, rn := 1:.N]
  using_annot = 'annot' %in% names(maf)

  if('problem' %in% names(maf) && maf[! is.na(problem), .N] > 0) {
    message("Found a \"problem\" column in input MAF. Variant records with problems will be excluded from output.")
    maf = maf[is.na(problem)]
  }
  
  if(maf[, .N] == 0) {
    stop('The MAF has no valid records.')
  }
  
  # In order to get access to getSeq() and clean_granges_for_cesa(), simplest
  # to require either ces.refset.hg38 or ces.refset.hg19
  if(identical(genome, 'hg38')) {
    if(! require('ces.refset.hg38')) {
      stop("To run with genome = \"hg38\", the reference data package ces.refset.hg38 must be installed.")
    }
    refset_env = ces.refset.hg38::ces.refset.hg38
  } else if(identical(genome, 'hg19')) {
    if(! require('ces.refset.hg19')) {
      msg = paste0("To run with genome = \"hg19\", the reference data package ces.refset.hg19 must be installed. ",
                   "Alternatively, use the chain_file argument of preload_maf() to convert the MAF to hg38.")
      stop(pretty_message(msg, emit = F))
    }
    refset_env = ces.refset.hg19::ces.refset.hg19
  } else {
    msg = paste0('Argument genome must be either "hg38" or "hg19". If your MAF data uses a different genome build, ',
                 'you can convert it using the chain_file argument of preload_maf().')
    stop(pretty_message(msg, emit = F))
  }
  
  for_ref_check = maf[variant_type == 'snv', .(chr = Chromosome, start = Start_Position, end = Start_Position, Reference_Allele)]
  if(for_ref_check[, .N] > 0) {
    if(for_ref_check[, .N] > 1000) {
      for_ref_check = for_ref_check[sample(1:.N, size = 1000, replace = F)]
    }
    stated_ref = for_ref_check$Reference_Allele
    
    ref_check_gr = tryCatch(clean_granges_for_cesa(gr = makeGRangesFromDataFrame(for_ref_check), 
                                          refset_env = refset_env,
                                          reduce_sort_strip = FALSE),
                            error = function(e) {
                              if(conditionMessage(e) %like% 'out-of-bound ranges') {
                                msg = paste0('Could not convert MAF to ', genome, ' GRanges. Most likely, ',
                                             'the genome build is incorrect. (Specify with genome = "hg19" or genome = "hg38".) ',
                                             'If you are unsure, run the MAF through preload_maf() ',
                                             'to investigate further.')
                                stop(pretty_message(msg, emit = F))
                              } else {
                                stop(conditionMessage(e))
                              }
                            }
    )
    actual_ref = as.character(BSgenome::getSeq(refset_env$genome, ref_check_gr))
    if(any(actual_ref != stated_ref)) {
      msg = paste0('A random check of MAF rows found that values in Reference_Allele do not always match the specified reference genome ',
                   '(', genome, '). Verify that the genome build is correct, or run preload_maf() to investigate reference mismatches.')
      stop(pretty_message(msg, emit = F))
    }
  }
  
  if('end' %in% names(maf)) {
    maf$end = NULL
  }
  maf[, end := Start_Position + nchar(Reference_Allele) - 1]
  
  # Complex variants (all of type "other") have comma-delimited variant_id. We'll get end position by checking all.
  get_end = function(variant_id) {
    all_id = strsplit(variant_id, ',')[[1]]
    all_id_pos = as.numeric(gsub('(^.*:)|(_.*$)', '', all_id))
    all_id_ref = gsub('(^.*_)|(>.*$)', '', all_id)
    last_pos = max(all_id_pos + nchar(all_id_ref) - 1)
    return(last_pos)
  }
  maf[variant_id %like% ',', end := sapply(variant_id, get_end)]
  
  maf_gr = GenomicRanges::makeGRangesFromDataFrame(
    maf[, .(seqnames = Chromosome, start = Start_Position, end, variant_id, variant_type)], 
    keep.extra.columns = TRUE)
  
  maf_gr = clean_granges_for_cesa(gr = maf_gr, refset_env = refset_env, reduce_sort_strip = FALSE)
  
  # Get gene annotations
  if(genome == 'hg38') {
    pathscore_cds_ranges = readRDS(system.file('PathScore/PathScore_CDS_ranges_hg38.rds', package = 'cancereffectsizeR'))
  } else {
    pathscore_cds_ranges = readRDS(system.file('PathScore/PathScore_CDS_ranges_hg19.rds', package = 'cancereffectsizeR'))
  }
  nearest = as.data.table(GenomicRanges::distanceToNearest(x = maf_gr, subject = pathscore_cds_ranges, select = "all"))[distance == 0]
  gene_info_from_gtf = as.data.table(mcols(pathscore_cds_ranges[nearest$subjectHits])[, c("hugo_symbol", "entrez_id")])
  nearest = cbind(nearest[, .(queryHits)], gene_info_from_gtf[, .(hugo_symbol, entrez_id)])
  
  # Sometimes an sbs overlaps the same CDS region twice due to multiple GRanges (commonly, multiple same-gene transcripts that overlap).
  # Uniquify to get one row per gene match.
  nearest = unique(nearest, by = c("queryHits", "hugo_symbol"))
  
  if(using_annot) {
    output = maf[, .(patient_id, rn, annot)]
    output = unique(output[nearest, .(hugo_symbol, entrez_id, patient_id, annot), on = c(rn = 'queryHits')])
  } else {
    output = maf[, .(patient_id, rn)]
    output = unique(output[nearest, .(hugo_symbol, entrez_id, patient_id), on = c(rn = 'queryHits')])
  }
  
  if(! is.null(file)) {
    fwrite(output, file = file, sep = '\t')
    pretty_message(paste0('Wrote PathScore file to ', file, '.'), black = F)
  }
  return(output)
}
