#' Posterior Simulation of Model Coefficients
#'
#' Forwards model input to posterior simulation functions for vector error
#' correction models.
#'
#' @param object an object of class 'gvecmodel', usually, a result of a
#' call to \code{\link{create_gvecmodel}} in combination with
#' \code{\link[bvartools]{add_priors}} and \code{\link[bvartools]{add_initial_values}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvecmodel'.
#' 
#' 
#' @export
#' @method add_posterior_coefficients gvecmodel
add_posterior_coefficients.gvecmodel <- function(object, ...){
  
  object[["submodels"]] <- lapply(object[["submodels"]], add_posterior_coefficients, ...)
  
  return(object)
}