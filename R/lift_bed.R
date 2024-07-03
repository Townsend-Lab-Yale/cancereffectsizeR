#' Convert BED intervals between genome builds
#'
#' Use this utility to convert BED intervals between genome coordinate systems using liftOver. Only
#' the chr/start/end fields of the input BED are used (strand is ignored). The output GRanges
#' will have no associated seqinfo.
#' 
#' A warning is given if the lifted intervals are less than 95\% of the size of the original
#' intervals. When the BED input represents sequencing target intervals, most of the input
#' intervals will usually lift successfully.
#' 
#' @param bed Pathname of a BED file, or a GRanges (typically loaded from a BED file with \code{rtracklayer::import.bed()}).
#' @param chain A UCSC-style chain file, or a Chain object (such as from \code{rtracklayer::import.chain()}).
#' @param outfile If not NULL, the returned GRanges will be saved to the specified path using \code{rtracklayer::export.chain()}.
#' @return GRanges representing lifted intervals from input \code{bed}.
#' @export
lift_bed = function(bed, chain, outfile = NULL) {
  if(rlang::is_scalar_character(bed)) {
    if(! file.exists(bed)) {
      stop("Specified BED file does not exist.")
    }
    bed = rtracklayer::import.bed(bed)
  } else if(! is(bed, 'GRanges')) {
    stop('Input bed should be the path to a BED file or a GRanges object.')
  }
  bed = BiocGenerics::unstrand(bed)
  
  if(rlang::is_scalar_character(chain)) {
    if(! file.exists(chain)) {
      stop('Specified chain file does not exist.')
    }
    chain = rtracklayer::import.chain(chain)
  } else if (! is(chain, 'Chain')) {
    stop('Input chain should be the path to a UCSC-style chain file or a Chain object.')
  }
  names(chain) = sub("^chr", "", names(chain))
  seqlevelsStyle(bed) = 'NCBI'
  lifted = sort(reduce(unlist(rtracklayer::liftOver(bed, chain))))
  seqlevelsStyle(lifted) = 'NCBI'
  prop = sum(width(lifted)) / sum(width(reduce(bed)))
  if(prop < .95) {
    percent = round(prop * 100, 1)
    msg = paste0('The lifted intervals cover ', percent, '% of the width of the original intervals.',
                 ' (For BED files representing sequencing target regions, typically most intervals',
                 ' will successfully lift between genome builds. If the percentage is very low,',
                 ' make sure you have the correct genome build for the input file.)')
    warning(pretty_message(msg, emit = F))
  }
  if(! is.null(outfile)) {
    if(! rlang::is_scalar_character(outfile)) {
      stop('outfile should be NULL or a scalar character indicating a valid file path.')
    }
    if(! outfile %like% '\\.bed(\\.gz)?$') {
      stop('outfile must end in .bed or .bed.gz')
    }
    if(! dir.exists(dirname(outfile))) {
      stop('Directory specified in outfile does not exist.')
    }
    if(file.exists(outfile)) {
      stop('Specified outfile already exists.')
    }
    rtracklayer::export.bed(lifted, outfile)
    message('Lifted BED intervals saved to ', outfile, '.')
  }
  if(is.null(outfile)) {
    return(lifted)
  } else {
    return(invisible(lifted))
  }
}
