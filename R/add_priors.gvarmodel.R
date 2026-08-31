#' Add Priors to Bayesian Models
#'
#' Adds prior specifications to a list of models by passing each element to
#' the respective method.
#'
#' @param object a list of models, usually, the output of a call to
#' \code{\link{create_gvarmodel}} in combination with \code{\link{add_submodels}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvarmodel'.
#' 
#' @export
#' @method add_priors gvarmodel
add_priors.gvarmodel <- function(object, ...){
  
  # This assumes that lower level objects are of class 'modellist', which
  # is handled by bvartools-functions.
  for (i in 1:length(object[["submodels"]])) {
    object[["submodels"]][[i]] <- add_priors(object[["submodels"]][[i]], ...)
  }
  
  return(object)
}