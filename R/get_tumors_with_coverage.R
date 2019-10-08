#' get_tumors_with_coverage
#' @param coverage list as formatted in the coverage slot of a CESAnalysis
#' @param locus a single locus (chr:start-end) in a GRanges object
#' @return a character vector of tumor identifiers
#' 

get_tumors_with_coverage = function(coverage, locus) {
  covered = coverage$samples_by_coverage[["exome"]]
  for (panel in names(coverage$samples_by_coverage)) {
    if (panel == "exome") {
      next
    }
    if (sum(GenomicRanges::countOverlaps(coverage$granges[[panel]], locus)) > 0) {
      covered = c(covered, coverage$samples_by_coverage[[panel]])
    }
  }
  return(covered)
}