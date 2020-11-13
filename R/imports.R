#' @import data.table
#' @import GenomeInfoDb
#' @import BSgenome
.datatable.aware=TRUE
.ces_ref_data = new.env()

# format a string the way R should automatically, then feed it to message()
pretty_message = function(msg, black = T) {
  msg = paste0(strwrap(msg), collapse = "\n")
  if (black) {
    msg = crayon::black(msg)
  }
  message(msg)
}


NULL
