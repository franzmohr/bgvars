#' Transform a VECX Sub-Model to a VARX Sub-Model
#'
#' Rewrites an estimated sub-model in error correction form as a sub-model in
#' levels.
#'
#' @param object an object of class 'vecxsubmodel'.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The transformation itself is
#' \code{\link[bvartools]{vec_to_var}}, which a 'vecxsubmodel' inherits through
#' its 'bvecmodel' class. What this method adds is the part of the
#' specification that describes a sub-model of a global model rather than a
#' Bayesian VAR --- how many of its regressors are its own variables, how many
#' are the weakly exogenous variables of the other units, and with how many lags
#' each of them enters. \code{\link{submodels_to_gvar}} reads exactly those
#' entries when it stacks the sub-models, so without them a converted sub-model
#' could not be combined into a global model.
#'
#' The coefficients of the result are ordered the way an estimated VARX
#' sub-model orders them: the lags of the endogenous variables, then the weakly
#' exogenous variables from lag zero, then the deterministic terms. A trend that
#' entered the error correction term restricted becomes an ordinary
#' deterministic term of the level representation.
#'
#' @return An object of class 'varxsubmodel'.
#'
#' @references
#'
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series
#' analysis} (2nd ed.). Berlin: Springer.
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
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      coint = list(v_i = 0, p_tau_i = 1),
#'                      sigma = list(df = 3, scale = 0.0001))
#' object <- add_initial_values(object)
#' object <- add_posterior_coefficients(object)
#'
#' # Rewrite one sub-model in levels
#' model <- vec_to_var(object[["submodels"]][["US"]][[1]])
#'
#' @export
#' @method vec_to_var vecxsubmodel
vec_to_var.vecxsubmodel <- function(object, ...) {

  specification <- object[["model"]]

  # Weakly exogenous variables of a GVEC sub-model that also uses global
  # variables are not built by create_vecxsubmodel in the first place. Rather
  # than let the position arithmetic below quietly assume they are not there,
  # say so.
  if (isTRUE(specification[["m_global"]] > 0)) {
    stop("Global variables in a GVEC sub-model are not supported yet.\n",
         "Feel free to send a feature request.")
  }

  # A 'bvecmodel' counts its weakly exogenous variables in 'm' and their lagged
  # differences in 's'. A GVEC sub-model does not: create_vecxsubmodel records
  # the whole differenced block in 'm' and leaves 's' at one. The two agree
  # whenever the weakly exogenous variables enter with a single lag, and part
  # company as soon as they do not, so the specification is translated here
  # rather than left to coincide.
  object[["model"]][["m"]] <- specification[["k_exogen"]]
  object[["model"]][["s"]] <- specification[["p_exogen"]]

  result <- NextMethod()

  # Everything the level representation says about itself is taken from the
  # transformation, which is what actually built the regressors. Only the
  # sub-model specification is added on top of it.
  model <- result[["model"]]

  model[["type"]] <- if (model[["k"]] == 1) "ARX" else "VARX"
  model[["k_endogen"]] <- model[["k"]]
  model[["p_endogen"]] <- model[["p"]]
  model[["k_exogen"]] <- length(model[["exogen"]])
  model[["p_exogen"]] <- model[["s"]]
  model[["m_global"]] <- 0L
  model[["s_global"]] <- 0L

  for (entry in c("global", "varsel", "iterations", "burnin", "algorithm")) {
    if (!is.null(specification[[entry]])) {
      model[[entry]] <- specification[[entry]]
    }
  }

  # The rank is a property of the error correction form and says nothing about
  # the model in levels, but it is what a sub-model of a GVEC model is chosen
  # by, so it is kept for reference.
  model[["rank"]] <- specification[["rank"]]

  result[["model"]] <- model

  class(result) <- c("varxsubmodel", "bvarmodel", "list")

  return(result)
}
