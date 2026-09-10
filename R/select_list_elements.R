#' Select List Elements
#' 
#' Limits the elements of a named list to the specified entities, preserving
#' all class information.
#' 
#' @param object a named list whose elements should be filtered.
#' @param submodel character vector of the elements in \code{object}, which
#' should remain in the output.
#' 
#' @export
select_list_elements <- function(object, submodel) {
  
  if (is.null(names(object))) {
    stop("Argument 'object' must be a named list.")
  }
  
  check_avail <- submodel %in% names(object)
  if (any(!check_avail)) {
    stop("The following submodels are not available in argument 'object': ",
         paste0(submodel[which(!check_avail)], collapse = ", "), ".")
  }
  
  orig_class <- class(object)
  object <- object[submodel]
  class(object) <- orig_class
  
  return(object)
}