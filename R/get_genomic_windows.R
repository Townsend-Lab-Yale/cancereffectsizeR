#' Split genome into same-width ranges.
#' 
#' @param window_size Integer size of windows (default 100,000 bp). This function will be slow if
#'   window_size is unexpectedly small.

#' @param exclude_ENCODE_blacklist TRUE/FALSE on filtering out windows in which at least 50% of
#'   bases are within ENCODE blacklist regions (default TRUE).
#' @param exclude_Y TRUE/FALSE on whether to exclude chrY from output ranges; default TRUE.
#' @param annotate_intervals Add some gene-based annotation columns (default FALSE).
get_genomic_windows = function(refset, window_size = 1e5, exclude_ENCODE_blacklist = TRUE,
                               exclude_Y = TRUE, annotate_intervals = FALSE) {
  
  data_dir = ifelse(is.environment(refset), refset$data_dir,
                    parse_refset_name(refset)$data_dir)
  genome_info = get_ref_data(data_dir_or_cesa = data_dir, 'genome_build_info')

#   Previously, an intervals argument. Don't need, currently.
#   intervals: If supplied, will only split up the input intervals, as opposed to the full
#     primary assembly reference ranges. Intervals should be supplied as a data.table with fields
#     chr/start/end.
  
  bsg = BSgenome::getBSgenome(genome_info$BSgenome)
  suppressWarnings({GenomeInfoDb::seqlevelsStyle(bsg) = genome_info$chromosome_style})
  chr_to_use = genome_info$supported_chr
  chr_lengths = seqlengths(bsg)[chr_to_use]
  intervals = data.table(chr = names(chr_lengths), start = 1, end = chr_lengths)
  
  if(! is.data.table(intervals)) {
    stop('intervals should be a data.table')
  }
  missing_col = setdiff(c('chr', 'start', 'end'), names(intervals))
  if(length(missing_col) > 0) {
    stop("Missing columns in intervals table: ", paste(missing_col, collapse = ', '), '.')
  }
  if(! is.character(intervals$chr)) {
    stop('chr column in intervals table should be type character')
  }
  if(! all(intervals$chr %in% genome_info$supported_chr)) {
    stop('The names of one or more chromosomes in the input intervals table are not among the refset\'s supported chromosomes.')
  }
  if(! rlang::is_integerish(intervals$start)) {
    stop('start field in intervals table should be numeric.')
  }
  if(! rlang::is_integerish(intervals$end)) {
    stop('end field in intervals table should be numeric.')
  }
  
  if(! rlang::is_bool(exclude_Y)) {
    stop('exclude_Y should be TRUE/FALSE.')
  }
  if(exclude_Y) {
    intervals = intervals[chr != 'Y']
  }
  
  if(! rlang::is_scalar_integerish(window_size, finite = TRUE) || window_size < 10 || window_size > 1e8) {
    stop('window_size should be scalar numeric on [10, 1e8]')
  }

  chr_ranges = rbindlist(mapply(function(chr, start, end) 
  { 
    dt = data.table(chr = chr, starting = seq(start, end - window_size, window_size))
    dt$ending = dt$start + window_size - 1
    return (dt)
  }, intervals$chr, intervals$start, intervals$end, SIMPLIFY = FALSE))
  
  if(! rlang::is_bool(exclude_ENCODE_blacklist)) {
    stop('exclude_ENCODE_blacklist should be TRUE/FALSE.')
  }
  if(exclude_ENCODE_blacklist) {
    if(! check_for_ref_data(data_dir, 'ENCODE_blacklist_v2')) {
      stop('exclude_ENCODE_blacklist = TRUE, but the reference data set does not include the required information.')
    }
    bl = as.data.table(get_ref_data(data_dir, 'ENCODE_blacklist_v2'))
    setnames(bl, 'seqnames', 'chr')
    setkey(bl, chr, start, end)
    
    # Removing ranges that are at least half overlapping blacklisted regions (~200Mb, cumulatively, for hg38)
    chr_ranges[, range_id := 1:.N]
    setkey(chr_ranges, chr, starting, ending)
    in_bl = foverlaps(chr_ranges, bl, nomatch = NULL, type = 'any')
    in_bl[, ol_size := pmin(end, ending) - pmax(start, starting) + 1]
    prop_in_bl = in_bl[, .(total_ol_prop = sum(ol_size)/sum(ending - starting + 1)), by = 'range_id']
    ranges_to_drop = prop_in_bl[total_ol_prop >= .5, range_id]
    chr_ranges = chr_ranges[! in_bl, on = c('chr', 'starting', 'ending')]
    chr_ranges[, range_id := NULL]
  }
  setnames(chr_ranges, c('starting', 'ending'), c('start', 'end'))

  if(! rlang::is_bool(annotate_intervals)) {
    stop('annotate_intervals should be TRUE/FALSE.')
  }
  
  if(annotate_intervals) {
    if(is.environment(refset) && ! is.null(refset$cancer_gene_coord)) {
      cg = copy(refset$cancer_gene_coord)
    } else if(! check_for_ref_data('gene_coord_cancer_anno')) {
      stop('annotate_intervals = TRUE, but the needed annotations are not available with this refset.')
    } else {
      cg = copy(get_ref_data(data_dir, 'gene_coord_cancer_anno'))
    }
    setkey(cg, chr, start, end)
    
    # Note some windows have no genes
    chr_ranges[, rn := 1:.N]
    gw_anno = foverlaps(chr_ranges, cg, nomatch = NA)
    gw_anno[, genes := paste0(gene, collapse = '|'), by = 'rn']
    gw_anno[, cancer_grp := fcase(cancer_anno %in% c('oncogene', 'TSG'), 1,
                                  cancer_anno == 'other', 2,
                                  default = 3)]
    gw_anno[, ol_size := pmin(i.end, end) - pmax(i.start, start)]
    final_anno = gw_anno[order(rn, cancer_anno, -ol_size)][, .SD[1], by = 'rn',
                                                           .SDcols = c('chr', 'i.start', 'i.end', 'gene', 'genes', 'cancer_anno')]
    setnames(final_anno, c('i.start', 'i.end', 'gene', 'genes', 'cancer_anno'), c('start', 'end', 'gene1', 'all_genes', 'gene1_anno'))
    final_anno[, rn := NULL]
    chr_ranges = final_anno
  }

  return(chr_ranges[])
}
