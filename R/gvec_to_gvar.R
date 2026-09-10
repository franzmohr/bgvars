#' Transform a Global VEC Model to a Global VAR Model
#'
#' Rewrites every estimated sub-model of a global model in levels, so that the
#' global model can be solved.
#'
#' @param object an object of class 'gvecmodel' with estimated sub-models.
#'
#' @details
#' \code{\link{submodels_to_gvar}} stacks sub-models in levels. A GVEC model
#' therefore has to be rewritten before it can be solved, which this function
#' does for all of its sub-models at once, using
#' \code{\link{vec_to_var.vecxsubmodel}} on each of them.
#'
#' The global data, the variable index and the weight matrices are left
#' untouched: they describe the units of the model and their weights, neither of
#' which the transformation changes.
#'
#' @return An object of class 'gvarmodel'.
#'
#' @examples
#'
#' # Load data
#' data("dees2007")
#' submodel_data <- dees2007[["submodel_data"]]
#' global_data <- dees2007[["global_data"]]
#'
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("US", "JP", "CA"))
#'
#' # Set up the model
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 1999:2001)
#' object <- add_submodels(object,
#'                         endogen = c("y", "Dp"), p_endogen = 1,
#'                         exogen = c("y", "Dp"), p_exogen = 1,
#'                         const = "unrestricted", trend = "restricted", r = 1,
#'                         iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' object <- align_model_obs(object)
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      coint = list(v_i = 0, p_tau_i = 1),
#'                      sigma = list(df = 3, scale = 0.0001))
#' object <- add_initial_values(object)
#' object <- add_posterior_coefficients(object)
#'
#' # Rewrite in levels and solve
#' object <- gvec_to_gvar(object)
#' gvar <- submodels_to_gvar(object)
#'
#' @export
gvec_to_gvar <- function(object) {

  if (!"gvecmodel" %in% class(object)) {
    stop("Argument 'object' must be of class 'gvecmodel'.")
  }

  if (is.null(object[["submodels"]])) {
    stop("Argument 'object' does not contain sub-models.")
  }

  for (submodel in names(object[["submodels"]])) {

    models <- object[["submodels"]][[submodel]]
    model_class <- class(models)

    for (i in seq_along(models)) {

      if (is.null(models[[i]][["posterior"]])) {
        stop("Sub-model ", submodel, " has not been estimated yet, so it ",
             "cannot be rewritten in levels.")
      }

      models[[i]] <- vec_to_var(models[[i]])
    }

    class(models) <- model_class
    object[["submodels"]][[submodel]] <- models
  }

  class(object) <- c("gvarmodel", "list")

  return(object)
}
