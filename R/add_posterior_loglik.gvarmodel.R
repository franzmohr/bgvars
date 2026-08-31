#' Add Log-Likelihood
#'
#' Calculates and adds posterior log-likelihoods to the sub-models of an object
#' of class 'gvarmodel'.
#'
#' @param object an object of class 'gvarmodel'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvarmodel'.
#' 
#' 
#' @export
#' @method add_posterior_loglik gvarmodel
add_posterior_loglik.gvarmodel <- function(object, ...){
  
  for (i in 1:length(object[["submodels"]])) {
    object[["submodels"]][[i]] <- add_posterior_loglik(object[["submodels"]][[i]], ...)
  }
  
  return(object)
}