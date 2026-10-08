#' Scale the Error Correction Terms of a GVEC Model
#'
#' Puts the series in the error correction term of every sub-model on a
#' comparable scale before posterior simulation, and back on the scale of the
#' data afterwards.
#'
#' @param object an object of class 'gvecmodel', usually, the output of a call to
#' \code{\link{create_gvecmodel}} in combination with \code{\link{add_submodels}}.
#' @param ... further arguments passed to
#' \code{\link[bvartools]{scale_error_correction}} and
#' \code{\link[bvartools]{rescale_error_correction}}, such as \code{scale} and
#' \code{centre}.
#'
#' @details
#' The series in the error correction term of a sub-model are on the scale of
#' the data, and those scales differ a great deal: a log level and an inflation
#' rate enter the same term, and so does a linear trend running into the
#' hundreds. The prior on the cointegration space is isotropic, so it treats all
#' of those series as if they were comparably scaled.
#'
#' \code{scale_error_correction} passes every sub-model to
#' \code{\link[bvartools]{scale_error_correction}}, which divides each stochastic
#' series by the standard deviation of its first differences and a trend by its
#' own standard deviation, and stores the factors in the attribute
#' \code{"scale"} of the element \code{data$train$w}.
#' \code{rescale_error_correction} passes every sub-model to
#' \code{\link[bvartools]{rescale_error_correction}}, which puts the series and
#' the draws of \eqn{\beta} back on the scale of the data and drops the
#' attribute.
#'
#' The transformation is exact: with \eqn{D} the diagonal matrix of scaling
#' factors, \eqn{\alpha \beta' w_{t-1}} is the same term as
#' \eqn{\alpha (D \beta)' (D^{-1} w_{t-1})}, so the model that is estimated is
#' the one that was specified. Only \eqn{\beta} is affected, \eqn{\alpha} is not.
#'
#' Starting values of \eqn{\beta} that a sub-model already carries follow the
#' same transformation -- for time varying cointegration both the path and the
#' state before the sample -- so the chain starts from the model it simulates.
#' Whether \code{\link[bayests]{add_initial_values}} is called before or after
#' scaling makes no difference to what comes back: rescaling returns the
#' starting values on the scale of the data either way.
#'
#' With \code{centre = TRUE}, which \code{scale_error_correction} passes on,
#' the sample mean of each stochastic series is subtracted from it as well, and
#' with \code{scale = FALSE, centre = TRUE} it is only centred. The unrestricted
#' constant of a sub-model takes up the means, so a sub-model needs one, and
#' \code{rescale_error_correction} shifts the draws of the constant back, period
#' by period for a time varying sub-model. For a sub-model whose cointegration
#' vectors vary over time, centring keeps a step of those vectors from acting
#' as a random walk intercept on log levels far from zero; see section 'Prior on
#' the cointegration space' of \code{\link{add_priors.vecxsubmodel}} and the
#' details of \code{\link[bvartools]{scale_error_correction}}. Centring alone
#' leaves the starting values of \eqn{\beta} as they are.
#'
#' Scaling and centring change what the priors say, so they are
#' modelling decisions rather than details of the sampler.
#' \code{\link[=add_posterior_coefficients.gvecmodel]{add_posterior_coefficients}}
#' therefore never scales on its own. A scaled model has to be rescaled before its
#' sub-models are transformed into level form or combined into a global model.
#'
#' @return An object of class 'gvecmodel'.
#'
#' @examples
#'
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
#'
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("EA", "US"))
#'
#' # Create empty model
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#'
#' # Add weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 1999:2001)
#'
#' # Create sub-models
#' object <- add_submodels(object,
#'                         p_endogen = 1, p_exogen = 2,
#'                         global = "poil",
#'                         s = 2,
#'                         r = 1,
#'                         const = "unrestricted", trend = "restricted",
#'                         iterations = 10,
#'                         burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' # Add priors and initial values
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      coint = list(v_i = 0, p_tau_i = 1),
#'                      sigma = list(df = 3, scale = 0.0001))
#' object <- add_initial_values(object)
#'
#' # Simulate on scaled error correction terms
#' object <- scale_error_correction(object)
#' object <- add_posterior_coefficients(object)
#' object <- rescale_error_correction(object)
#'
#' @export
#' @method scale_error_correction gvecmodel
scale_error_correction.gvecmodel <- function(object, ...) {

  .map_vecxsubmodels(object, function(model) {
    model <- scale_error_correction(model, ...)
    # Centring alone leaves beta in the units of the series, so only scaling
    # carries factors to transform the starting values by.
    factors <- attr(model[["data"]][["train"]][["w"]], "scale")
    if (is.null(factors)) {
      return(model)
    }
    .transform_initial_beta(model, factors)
  })
}

#' @rdname scale_error_correction.gvecmodel
#' @export
#' @method rescale_error_correction gvecmodel
rescale_error_correction.gvecmodel <- function(object, ...) {

  .map_vecxsubmodels(object, function(model) {
    # Read before the call, which drops the attribute -- and which refuses a
    # sub-model that was neither scaled nor centred. A centred sub-model that
    # was not scaled has no factors, and its starting values of beta stay.
    factors <- attr(model[["data"]][["train"]][["w"]], "scale")
    model <- rescale_error_correction(model, ...)
    if (is.null(factors)) {
      return(model)
    }
    .transform_initial_beta(model, 1 / factors)
  })
}

# Apply a function to every model of every sub-model of a 'gvecmodel', keeping
# the class of the list each sub-model's models are kept in.
.map_vecxsubmodels <- function(object, fun) {

  object[["submodels"]] <- lapply(object[["submodels"]], function(models) {
    model_class <- class(models)
    models <- lapply(models, fun)
    class(models) <- model_class
    models
  })

  return(object)
}

# Multiply the starting values of beta by the diagonal matrix of 'factors'.
#
# add_initial_values() stores them stacked -- a column of k_ect * r elements, or
# of k_ect * r * t for the path of time varying cointegration vectors, with the
# state before the sample as 'beta_init' -- so they are reshaped to k_ect rows
# before the factors are applied, and written back into the object they came
# from, which keeps whatever dimensions it had. diag() is given an explicit
# size, because for a single factor it would otherwise read its argument as the
# size of an identity matrix rather than as the diagonal to build.
.transform_initial_beta <- function(model, factors) {

  scaling <- diag(factors, nrow = length(factors))
  for (name in c("beta", "beta_init")) {
    value <- model[["initial"]][[name]]
    if (!is.null(value)) {
      value[] <- as.numeric(scaling %*% matrix(value, length(factors)))
      model[["initial"]][[name]] <- value
    }
  }

  return(model)
}
