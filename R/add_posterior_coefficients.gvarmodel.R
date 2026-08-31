#' Posterior Simulation of Model Coefficients
#'
#' Forwards model input to posterior simulation functions for vector autoregressive models.
#'
#' @param object an object of class 'gvarmodel', usually, a result of a
#' call to \code{\link{create_gvarmodel}} in combination with
#' \code{\link[bvartools]{add_priors}} and \code{\link[bvartools]{add_initial_values}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvarmodel'.
#' 
#' 
#' @export
#' @method add_posterior_coefficients gvarmodel
add_posterior_coefficients.gvarmodel <- function(object, ...){
  
  object[["submodels"]] <- lapply(object[["submodels"]], add_posterior_coefficients, ...)
  
  return(object)
}