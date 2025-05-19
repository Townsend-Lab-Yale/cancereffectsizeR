#' Split genome into same-width ranges.
#' 
#' ## TO-DO: change to excluding windows that are >50% blacklisted, instead of the weird criteria here.
#' @param window_size Integer size of windows. Will be a little slow if windows are very small.
#' @param ranges If supplied, will only split up the input GRanges, as opposed to the full primary assembly reference ranges.
#' @param exclude_ENCODE_blacklist TRUE/FALSE on filtering out windows that overlap large (> .5 * window_size)
#' regions flagged in ENCODE blacklist v2 (default TRUE).
#' @param exclude_non_PAR TRUE/FALSE on removing windows that overlap non-pseudoautosomal regions of chrX/chrY (default TRUE).
get_genomic_windows = function(refset, window_size = 1e5, ranges = NULL, exclude_ENCODE_blacklist = TRUE,
                               exclude_non_PAR = TRUE) {
  
  data_dir = ifelse(is.environment(refset), refset$data_dir,
                    parse_refset_name(refset)$data_dir)
  
  ## To-do: add support for NULL ranges
  # genome_info = get_ref_data(data_dir_or_cesa = data_dir, 'genome_build_info')
  # bsg = BSgenome::getBSgenome(genome_info$BSgenome)
  # suppressWarnings({GenomeInfoDb::seqlevelsStyle(bsg) = genome_info$chromosome_style})
  if(is.null(ranges)) {
    stop('ranges is NULL. Support for this is coming.')
  }
  
  
  if(exclude_ENCODE_blacklist) {
    bl = as.data.table(get_ref_data(data_dir, 'ENCODE_blacklist_v2'))
    setnames(bl, 'seqnames', 'chr')
  }
  
  # To-do: remove hardcoding
  
  # confirmed via getSeq for hg38
  par1_end = 2781479
  par2_start = 155701383
  
  # Ensembl blacklist regions
  
  # We'll ignore small regions
  bl = bl[width > window_size / 2]
  setkey(bl, chr, start, end)
  
  
  if(! rlang::is_scalar_integerish(window_size, finite = TRUE) || window_size < 10 || window_size > 1e8) {
    stop('window_size should be scalar numeric on [10, 1e8]')
  }
  

  if(! is.data.table(ranges)) {
    stop('ranges should be a data.table')
  }
  chr_ranges = rbindlist(mapply(function(chr, start, end) 
  { 
    dt = data.table(chr = chr, starting = seq(start, end - window_size, window_size))
    dt$ending = dt$start + window_size - 1
    return (dt)
  }, ranges$chr, ranges$p_start, ranges$q_end, SIMPLIFY = FALSE))
  
  # For now, leave out chrX non-PAR regions
  chr_ranges = chr_ranges[!(chr == 'X' & (starting < par1_end | ending > par2_start))]
  
  # Also removing Ensembl blacklist regions
  setkey(chr_ranges, chr, starting, ending)
  in_bl = foverlaps(chr_ranges, bl, nomatch = NULL)
  chr_ranges = chr_ranges[! in_bl, on = c('chr', 'starting', 'ending')]
  return(chr_ranges)
}
