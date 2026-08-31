#' Add Log-Likelihood
#'
#' Calculates and adds posterior log-likelihoods to the sub-models of an object
#' of class 'gvecmodel'.
#'
#' @param object an object of class 'gvecmodel'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return A list of class 'gvecmodel'.
#' 
#' 
#' @export
#' @method add_posterior_loglik gvecmodel
add_posterior_loglik.gvecmodel <- function(object, ...){
  
  for (i in 1:length(object[["submodels"]])) {
    object[["submodels"]][[i]] <- add_posterior_loglik(object[["submodels"]][[i]], ...)
  }
  
  return(object)
}