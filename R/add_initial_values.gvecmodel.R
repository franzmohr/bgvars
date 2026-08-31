#' Add Initial Values of an MCMC Chain
#'
#' Adds initial values to a list of models by passing each element to
#' the respective method.
#'
#' @param object an object of class 'gvecmodel', usually, the output of a call to
#' \code{\link{create_gvecmodel}} in combination with \code{\link{add_submodels}}
#' and \code{\link[bvartools]{add_initial_values}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvecmodel'.
#' 
#' 
#' @export
#' @method add_initial_values gvecmodel
add_initial_values.gvecmodel <- function(object, ...){
  
  # This assumes that lower level objects are of class 'modellist', which
  # is handled by bvartools-functions.
  for (i in 1:length(object[["submodels"]])) {
    object[["submodels"]][[i]] <- add_initial_values(object[["submodels"]][[i]], ...)
  }
  
  return(object)
}