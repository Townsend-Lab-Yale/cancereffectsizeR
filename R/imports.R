#' @import data.table
#' @import GenomeInfoDb
#' @import BSgenome
.datatable.aware=TRUE
.ces_ref_data = new.env()
.official_ref_sets = list(ces.refset.hg19 = as.package_version("0.0.0.9000"))

# format a string the way R should automatically, then feed it to message()
pretty_message = function(msg, emit = T, black = emit) {
  msg = paste0(strwrap(msg), collapse = "\n")
  if (black) {
    msg = crayon::black(msg)
  }
  if (emit) {
    message(msg)
  } else {
    return(msg)
  }
}


NULL
