#' Posterior Simulation of Model Coefficients
#'
#' Forwards model input to posterior simulation functions for vector error
#' correction models.
#'
#' @param object an object of class 'gvecmodel', usually, a result of a
#' call to \code{\link{create_gvecmodel}} in combination with
#' \code{\link[bayests]{add_priors}} and \code{\link[bayests]{add_initial_values}}.
#' @param ... further arguments passed to or from other methods.
#' @param cores the number of worker processes the models are simulated on.
#' With more than one, all models of all sub-models are shared out over one
#' cluster, as \code{\link[bayests]{add_posterior_coefficients}} does for a
#' 'modellist'. The draws do not depend on the number of workers, since every
#' model draws with its own seed (see \code{\link[=add_seed.gvarmodel]{add_seed}}).
#'
#' @details
#' Every sub-model is simulated on the error correction term it carries. The
#' function does not change the scale of that term.
#'
#' The series in the error correction term of a sub-model are on the scale of
#' the data, and those scales differ a great deal, while the prior on the
#' cointegration space is isotropic. To simulate on a comparable scale, call
#' \code{\link[=scale_error_correction.gvecmodel]{scale_error_correction}} before
#' this function and
#' \code{\link[=scale_error_correction.gvecmodel]{rescale_error_correction}}
#' after it:
#' \preformatted{object <- scale_error_correction(object)
#' object <- add_posterior_coefficients(object)
#' object <- rescale_error_correction(object)}
#' Scaling changes what the prior on the cointegration space says, so it is left
#' to the user as a step of its own rather than done implicitly.
#'
#' @return A list of class 'gvecmodel'.
#'
#'
#' @export
#' @method add_posterior_coefficients gvecmodel
add_posterior_coefficients.gvecmodel <- function(object, ..., cores = 1){

  # The argument used to switch an implicit scaling on and off. Passed on
  # through '...' it would be ignored by the samplers without a word, and a
  # script written for the old default would simulate on unscaled terms without
  # noticing, so it is refused.
  if ("scale_ect" %in% names(list(...))) {
    stop("Argument 'scale_ect' is no longer supported. Call 'scale_error_correction' ",
         "before and 'rescale_error_correction' after 'add_posterior_coefficients' ",
         "to simulate on scaled error correction terms.")
  }

  .apply_to_submodel_models(object, add_posterior_coefficients, ..., cores = cores)
}
