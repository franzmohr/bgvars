#' Posterior Simulation of Model Coefficients
#'
#' Forwards model input to posterior simulation functions for vector error
#' correction models.
#'
#' @param object an object of class 'gvecmodel', usually, a result of a
#' call to \code{\link{create_gvecmodel}} in combination with
#' \code{\link[bvartools]{add_priors}} and \code{\link[bvartools]{add_initial_values}}.
#' @param scale_ect logical. Should the series in the error correction term be
#' put on a comparable scale for the posterior simulation? Defaults to
#' \code{TRUE}. See 'Details'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The series in the error correction term of a sub-model are on the scale of
#' the data, and those scales differ a great deal: a log level and an inflation
#' rate of the same unit enter the same term, and so does a linear trend running
#' into the hundreds. The prior on the cointegration space is isotropic --
#' \code{\link[bvartools]{add_priors}} builds it from \code{coint$p_tau_i} as a
#' multiple of the identity matrix -- so it treats all of those series as if they
#' were comparably scaled.
#'
#' With \code{scale_ect = TRUE} each sub-model is therefore simulated on a
#' scaled error correction term, using
#' \code{\link[bvartools]{scale_error_correction}}, and the draws are put back
#' on the scale of the data afterwards with
#' \code{\link[bvartools]{rescale_error_correction}}. The transformation is
#' exact: with \eqn{D} the diagonal matrix of scaling factors,
#' \eqn{lpha eta' w_{t-1}} is unchanged by writing it as
#' \eqn{lpha (D eta)' (D^{-1} w_{t-1})}, so the model that is estimated is
#' the same one. Only the draws of \eqn{eta} are affected by the
#' transformation, those of \eqn{lpha} are not.
#'
#' What the function returns is on the scale of the data either way. The initial
#' values the object carries are left untouched, so that the reduced rank
#' maximum likelihood estimates \code{\link[bvartools]{add_initial_values}}
#' produced remain comparable with published cointegrating vectors.
#'
#' @return A list of class 'gvecmodel'.
#'
#'
#' @export
#' @method add_posterior_coefficients gvecmodel
add_posterior_coefficients.gvecmodel <- function(object, scale_ect = TRUE, ...){

  if (!is.logical(scale_ect) || length(scale_ect) != 1) {
    stop("Argument 'scale_ect' must be a single logical value.")
  }

  draw <- if (scale_ect) .draw_scaled_vecxsubmodel else add_posterior_coefficients

  object[["submodels"]] <- lapply(object[["submodels"]], function(models) {
    model_class <- class(models)
    models <- lapply(models, draw, ...)
    class(models) <- model_class
    models
  })

  return(object)
}