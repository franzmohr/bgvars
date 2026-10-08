#' Posterior Simulation of Model Coefficients
#'
#' Forwards model input to posterior simulation functions for vector autoregressive models.
#'
#' @param object an object of class 'gvarmodel', usually, a result of a
#' call to \code{\link{create_gvarmodel}} in combination with
#' \code{\link[bayests]{add_priors}} and \code{\link[bayests]{add_initial_values}}.
#' @param ... further arguments passed to or from other methods.
#' @param cores the number of worker processes the models are simulated on.
#' With more than one, all models of all sub-models are shared out over one
#' cluster, as \code{\link[bayests]{add_posterior_coefficients}} does for a
#' 'modellist'. The draws do not depend on the number of workers, since every
#' model draws with its own seed (see \code{\link[=add_seed.gvarmodel]{add_seed}}).
#' 
#' @return A list of class 'gvarmodel'.
#' 
#' 
#' @export
#' @method add_posterior_coefficients gvarmodel
add_posterior_coefficients.gvarmodel <- function(object, ..., cores = 1){
  .apply_to_submodel_models(object, add_posterior_coefficients, ..., cores = cores)
}
