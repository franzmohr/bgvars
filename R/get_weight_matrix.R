#' Get Weight Matrix
#'
#' Obtains the weight matrix of a sub-model, which links the vector of all
#' endogenous variables of the global model to the endogenous and weakly
#' exogenous variables of that sub-model.
#'
#' @param object an object of class 'gvarmodel' or 'gvecmodel' containing
#' weight matrices, that is, the output of a call to
#' \code{\link{add_weight_matrices}}.
#' @param submodel character of the sub-model, whose weight matrix should be
#' obtained.
#' @param period integer of the period, for which the weight matrix should be
#' obtained. Defaults to 1. Only relevant for time varying weights.
#'
#' @details
#' \code{\link{add_weight_matrices}} stores the weight matrices of a sub-model
#' as a single matrix, in which the matrices of the individual periods are
#' stacked on top of each other. This function cuts out the matrix
#' \eqn{W_{i}} of one period and labels it, so that
#' \deqn{z_{it} = W_{i} y_{t}}
#' with \eqn{z_{it} = (x_{it}', x_{it}^{*'})'} the vector of the endogenous and
#' weakly exogenous variables of sub-model \eqn{i} and \eqn{y_t} the vector of
#' the endogenous variables of the global model.
#'
#' The rows of the result are named after the variables of the sub-model, its
#' endogenous variables first, followed by its weakly exogenous variables. Note
#' that these names are not unique, since a variable usually enters a sub-model
#' both as an endogenous and as a weakly exogenous variable. The columns are
#' named after the variables of the global model.
#'
#' @return A matrix.
#'
#' @examples
#'
#' # Load data
#' data("gvar2019")
#' global_data <- gvar2019[["global_data"]]
#' submodel_data <- gvar2019[["submodel_data"]]
#'
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#'
#' # Create global model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#'
#' # Generate weight matrices
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 3)
#'
#' # Obtain the weight matrix of the Austrian sub-model
#' w <- get_weight_matrix(object, "AT")
#'
#' @export
get_weight_matrix <- function(object, submodel, period = 1) {

  if (is.null(object[["weights"]])) {
    stop("Argument 'object' does not contain weight matrices. See 'add_weight_matrices'.")
  }

  if (length(submodel) != 1) {
    stop("Argument 'submodel' may only contain one element.")
  }

  if (!submodel %in% names(object[["weights"]])) {
    stop("No weight matrix available for sub-model ", submodel, ".")
  }

  index <- object[["global"]][["index"]]
  vars_endogen <- index[which(index[, "submodel"] == submodel), "variable"]
  vars_exogen <- unique(index[index[, "submodel"] != submodel, "variable"])
  n_vars <- length(vars_endogen) + length(vars_exogen)

  weights <- object[["weights"]][[submodel]]

  # The matrices of the individual periods are stacked on top of each other.
  n_periods <- nrow(weights) / n_vars
  if (length(period) != 1 || period < 1 || period > n_periods) {
    stop("Argument 'period' must be a single integer between 1 and ",
         n_periods, ".")
  }

  temp <- weights[(period - 1) * n_vars + 1:n_vars, , drop = FALSE]
  dimnames(temp) <- list(c(vars_endogen, vars_exogen), index[, "index"])

  return(temp)
}

# The weight matrix of a sub-model, arranged the way the estimated coefficients
# of that sub-model are.
#
# get_weight_matrix() returns the rows in the order of the global index, which
# is the order of the data. A sub-model, however, orders its variables the way
# they were named in add_submodels(), and its coefficient matrices follow that
# order. The rows are therefore reordered here to match, and restricted to the
# variables the sub-model actually uses.
#
# 'keep' are the columns of the global variable vector that survive into the
# global model. Variables that no sub-model uses as an endogenous variable are
# not part of that vector, so their columns are dropped -- which is only sound
# if they carry no weight for the rows that remain. That is checked rather than
# assumed.
.submodel_weight_matrix <- function(object, submodel, period, keep) {

  model <- object[["submodels"]][[submodel]][[1]][["model"]]

  index <- object[["global"]][["index"]]
  vars_endogen <- index[which(index[, "submodel"] == submodel), "variable"]
  vars_exogen <- unique(index[index[, "submodel"] != submodel, "variable"])

  w <- get_weight_matrix(object, submodel, period = period)

  rows_endogen <- match(model[["endogen"]], vars_endogen)
  if (anyNA(rows_endogen)) {
    stop("For sub-model ", submodel, " the endogenous variable(s) ",
         paste0(model[["endogen"]][is.na(rows_endogen)], collapse = ", "),
         " are not available in the global variable index.")
  }

  # Weakly exogenous variables carry the suffix that create_varxsubmodel and
  # create_vecxsubmodel append to them.
  exogen <- sub("\\.s$", "", model[["exogen"]])
  rows_exogen <- match(exogen, vars_exogen)
  if (anyNA(rows_exogen)) {
    stop("For sub-model ", submodel, " the weakly exogenous variable(s) ",
         paste0(exogen[is.na(rows_exogen)], collapse = ", "),
         " are not available in the global variable index.")
  }
  rows_exogen <- length(vars_endogen) + rows_exogen

  rows <- c(rows_endogen, rows_exogen)

  dropped <- setdiff(seq_len(ncol(w)), keep)
  if (length(dropped) > 0) {
    if (any(w[rows, dropped, drop = FALSE] != 0)) {
      stop("Sub-model ", submodel, " puts weight on variables that are not ",
           "endogenous to any sub-model of the global model, so the global ",
           "model cannot be closed.\n",
           "Affected variables: ",
           paste0(index[dropped, "index"][
             apply(w[rows, dropped, drop = FALSE] != 0, 2, any)],
             collapse = ", "), ".")
    }
  }

  return(w[rows, keep, drop = FALSE])
}
