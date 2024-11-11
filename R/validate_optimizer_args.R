#' Check custom optimizer arguments
#'
#' @param optimizer_args List of arguments/values to pass to the optimizer.
#' @keywords internal
validate_optimizer_args = function(optimizer_args) {
  if (! is(optimizer_args, "list") || uniqueN(names(optimizer_args)) != length(optimizer_args)) {
    stop("optimizer_args should a named list of arguments to pass.")
  }
  
  reserved_args = c('minuslogl', 'start', 'vecpar')
  if(any(reserved_args %in% names(optimizer_args))) {
    msg = paste0('Optimizer arguments start, vecpar, and minuslogl cannot be changed here. ', 
                 'If you are using a custom model, your likelihood function can declare these ',
                 'values directly (see docs).')
    stop(pretty_message(msg, emit = F))
  }
  if('control' %in% names(optimizer_args) && ! is(optimizer_args$control, 'list')) {
    stop('Optimizer argument \"control\" should be a list.')
  }
}
