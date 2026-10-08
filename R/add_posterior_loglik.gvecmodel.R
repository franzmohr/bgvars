#' Add Log-Likelihood
#'
#' Calculates and adds posterior log-likelihoods to the sub-models of an object
#' of class 'gvecmodel'.
#'
#' @param object an object of class 'gvecmodel'.
#' @param ... further arguments passed to or from other methods.
#' @param cores the number of worker processes the log-likelihoods are
#' computed on. With more than one, all models of all sub-models are shared out
#' over one cluster, as \code{\link[bayests]{add_posterior_loglik}} does for
#' a 'modellist'.
#' 
#' @return A list of class 'gvecmodel'.
#' 
#' 
#' @export
#' @method add_posterior_loglik gvecmodel
add_posterior_loglik.gvecmodel <- function(object, ..., cores = 1){
  
  .apply_to_submodel_models(object, add_posterior_loglik, ..., cores = cores)
}
