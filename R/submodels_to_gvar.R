#' Generate Global Model from Sub-Models
#'
#' Combines the estimated sub-models of a global model into the solved global
#' model.
#'
#' @param object an object of class 'gvarmodel' or 'gvecmodel' containing
#' exactly one estimated model per sub-model, or one of class 'gvarexport', the
#' handle of an export. Use the \code{submodels} argument of
#' \code{\link{read_gvar_from_folder}}, or of this function for an export, to
#' reduce a model with several candidates per sub-model to one.
#' @param submodels for an export; a data frame with the columns
#' \code{submodel} and \code{position}, naming which model to take per
#' sub-model, as \code{\link{read_gvar_from_folder}} takes it. Not needed where
#' the export holds one model per sub-model.
#' @param chunk for an export; how many draws are read and stacked at a time.
#' Defaults to 100.
#' @param ... arguments passed to the method.
#' @param period integer of the period of the global data, at which the model
#' should be solved. Defaults to the last period of the estimation sample. It
#' selects the weight matrices, if they vary over time, and the coefficients and
#' error covariance matrices of sub-models with time varying parameters or
#' stochastic volatility, which must then lie within the estimation sample.
#'
#' @details
#' Each sub-model \eqn{i} is a VARX* model in its endogenous variables
#' \eqn{x_{it}} and its weakly exogenous variables \eqn{x^{*}_{it}}, which
#' are weighted averages of the endogenous variables of the remaining
#' sub-models. Written in terms of the vector \eqn{y_t} of the endogenous
#' variables of all sub-models, the contemporaneous terms of sub-model \eqn{i}
#' become \eqn{(A_{0i}, -A^{*}_{0i}) W_i y_t}, and stacking these blocks over
#' all sub-models gives the global matrix \eqn{G}. The lagged terms are
#' stacked the same way, which turns the sub-models into a single VAR in
#' \eqn{y_t}.
#'
#' The sub-models are stacked in levels. A sub-model in error correction form
#' is therefore rewritten as the level VARX model it implies first, with
#' \code{\link{vec_to_var.vecxsubmodel}}, so a GVEC model is solved in the same
#' single step as a GVAR model. The same holds for a model whose sub-models
#' were chosen from VARX and VECX candidates alike and so come in both forms.
#' The cointegration rank of a sub-model describes its error correction form
#' only and plays no further role. A time varying sub-model in error correction
#' form is reduced to the draws of \code{period} before it is rewritten, so only
#' the draws of that period are transformed.
#'
#' The model is returned in its reduced form, that is, premultiplied by
#' \eqn{G^{-1}}, which is the form the methods of \pkg{bvartools} work on.
#' \eqn{G} is a dense matrix and can therefore not be passed on as the
#' structural form of a 'bvarmodel', which stores the free elements of a unit
#' lower triangular matrix. Its draws are returned as element \code{g}
#' instead.
#'
#' The sub-models may use different deterministic terms. The global model uses
#' every term that any of them uses, in the order of the deterministic series of
#' the global model, and a sub-model has coefficients of zero on the terms it
#' does not use. This typically arises in a GVEC model, when a sub-model with a
#' restricted trend is estimated at cointegration rank zero: without an error
#' correction term the trend drops out of its level representation (see
#' \code{\link{vec_to_var.vecxsubmodel}}), while the sub-models of higher rank
#' keep it.
#'
#' The error covariance matrix of the structural form is block diagonal, one
#' block per sub-model. The reduced form errors \eqn{G^{-1} u_t} are
#' correlated across sub-models accordingly.
#'
#' Sub-models with time varying parameters or stochastic volatility hold one set
#' of coefficients and one error covariance matrix per period of the estimation
#' sample. The global model is solved from the draws of the period given in
#' \code{period}, so the result describes the global model at that date, and a
#' model for every period of interest is obtained by calling the function once
#' per period. Sub-models with constant coefficients and a constant error
#' covariance matrix contribute the same draws in every period.
#'
#' The reduced form has constant coefficients, a full error precision and no
#' structural form whatever the sub-models had, so the returned object names
#' \code{VarNormalWishart} as the algorithm its draws came from. That is what
#' the methods of \pkg{bvartools} that work on draws dispatch on, forecasting
#' included.
#'
#' @section Solving from an export:
#'
#' Stacking the sub-models is the one step of the package that needs every
#' sub-model at once, and it needs them one draw at a time: a draw of the global
#' model is built from that draw of every sub-model and from nothing else. Given
#' the handle of an export, the function reads the draws in pieces of
#' \code{chunk}, stacks each piece and puts the pieces end to end, so what it
#' holds of the sub-models is one piece rather than every chain. For a global
#' model of time varying sub-models that is the difference between a few
#' gigabytes and tens of them.
#'
#' The result is the one the same models in memory give, draw for draw. The
#' global model's own draws are held whole, which is what every method that
#' works on them expects. They are the smaller half -- one square matrix per lag
#' and draw, against a coefficient path per period and draw in every time
#' varying sub-model -- but for a large global model they are still worth
#' putting away: written with \code{\link[bvartools]{write_to_hdf5}} and opened
#' with \code{\link[bvartools]{open_model}}, the impulse responses and variance
#' decompositions of \pkg{bvartools} read them in pieces as well.
#'
#' \preformatted{
#' solved <- submodels_to_gvar(export)
#' write_to_hdf5(solved, filename = file)
#' rm(solved)
#'
#' stored <- open_model(file)
#' irf(stored, impulse = "US_r", response = "DE_y", n_ahead = 20, type = "gir")
#' }
#'
#' @return An object of class 'bvarmodel', which can be analysed with the
#' functionality of \pkg{bvartools} -- impulse responses, forecast error
#' variance decompositions and forecasts. A forecast needs the future path of
#' the global variables, which \code{\link[bvartools]{add_forecast_input}}
#' takes as its \code{exogen} argument; \code{data$original$exogen} of the
#' returned object holds the series to continue.
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
#' # Set up the model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 3)
#' object <- add_submodels(object,
#'                         p_endogen = 1,
#'                         exogen = c("y", "Dp", "eq", "r", "lr"),
#'                         p_exogen = 2,
#'                         global = "poil", s = 2,
#'                         error = "wishart",
#'                         iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' object <- align_model_obs(object)
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      sigma = list(df = 3, scale = 0.0001))
#' object <- add_initial_values(object)
#' object <- add_posterior_coefficients(object)
#'
#' # Solve the global model
#' gvar <- submodels_to_gvar(object)
#'
#' # Forecast it, holding the global variable at its last observed value
#' poil <- gvar[["data"]][["original"]][["exogen"]]
#' tsp_poil <- stats::tsp(poil)
#' future <- stats::ts(matrix(rep(poil[nrow(poil), ], each = 5), 5, ncol(poil)),
#'                     start = tsp_poil[2], frequency = tsp_poil[3])
#' dimnames(future)[[2]] <- dimnames(poil)[[2]]
#'
#' gvar <- add_forecast_input(gvar, n_ahead = 4, exogen = future)
#' gvar <- add_posterior_forecasts(gvar)
#' fcst <- predict(gvar, n_ahead = 4)
#'
#' # A GVEC model is solved the same way
#' data("dees2007")
#' submodel_data <- select_list_elements(dees2007[["submodel_data"]],
#'                                       c("US", "JP", "CA"))
#'
#' object <- create_gvecmodel(submodel_data = submodel_data,
#'                            global_data = dees2007[["global_data"]])
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 1999:2001)
#' object <- add_submodels(object,
#'                         endogen = c("y", "Dp"), p_endogen = 1,
#'                         exogen = c("y", "Dp"), p_exogen = 2,
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
#' gvec <- submodels_to_gvar(object)
#'
#' @export
submodels_to_gvar <- function(object, ...) {
  UseMethod("submodels_to_gvar")
}

#' @rdname submodels_to_gvar
#' @export
#' @method submodels_to_gvar gvarmodel
submodels_to_gvar.gvarmodel <- function(object, period = NULL, ...) {
  .combine_submodels(object, period = period)
}

#' @rdname submodels_to_gvar
#' @export
#' @method submodels_to_gvar gvecmodel
submodels_to_gvar.gvecmodel <- function(object, period = NULL, ...) {
  .combine_submodels(object, period = period)
}

#' @rdname submodels_to_gvar
#' @export
#' @method submodels_to_gvar gvarexport
submodels_to_gvar.gvarexport <- function(object, period = NULL, submodels = NULL,
                                         chunk = 100, ...) {
  .solve_export(object, period = period, submodels = submodels, chunk = chunk)
}

# What to say when a sub-model has a posterior but no draws.
#
# The discounted models write a closed form -- a mean, a scale and a covariance
# per period -- and deliberately no /posterior/a/coeffs: joining one i.i.d. draw
# per period would look like a sampled coefficient path and is not one, the
# smoothed posterior being dependent across periods. A global model is built by
# stacking one draw of every sub-model at a time, so there is nothing here to
# stack.
#
# Drawing from the closed form is what a host that wants draws does, and it is
# BayesTS's to do: the arithmetic is in VarTvpDiscountEstimator::draw_period(),
# and reimplementing it here would put the same posterior in two places.
.refuse_discount_submodel <- function(algorithm) {

  stop("Sub-models estimated with ", algorithm, " carry a closed-form posterior ",
       "rather than draws -- a mean, a scale and a covariance per period, and no ",
       "/posterior/a/coeffs -- so there is nothing to stack into a global model. ",
       "They are for comparing specifications by the exact log marginal ",
       "likelihood their files carry; use a sampling algorithm for the model the ",
       "impulse responses are taken from.")
}

# The algorithm of a sub-model that is one of the discounted models, or NULL.
.discount_submodel_algorithm <- function(specifications) {

  algorithms <- unique(stats::na.omit(specifications))
  discounted <- algorithms[vapply(algorithms, bvartools::is_discount_model,
                                  logical(1))]
  if (length(discounted) == 0) {
    return(NULL)
  }

  paste(discounted, collapse = " and ")
}

# The global model from the sub-models of an object that holds their draws.
.combine_submodels <- function(object, period = NULL) {

  if (is.null(object[["submodels"]])) {
    stop("Argument 'object' does not contain sub-models.")
  }

  algorithm <- .discount_submodel_algorithm(unlist(lapply(
    object[["submodels"]], function(models) {
      lapply(models, function(model) model[["model"]][["algorithm"]])
    })))
  if (!is.null(algorithm)) {
    .refuse_discount_submodel(algorithm)
  }

  index <- object[["global"]][["index"]]
  submodels <- unique(index[, "submodel"])

  # Number of models per sub-model
  n_models_i <- unlist(lapply(object[["submodels"]], function(x){length(x)}))
  if (any(n_models_i > 1)) {
    pos <- which(n_models_i > 1)
    msg_names <- paste0(names(n_models_i)[pos], collapse = ", ")
    stop("Argument 'object' may only contain one model per sub-model.\n",
         "Affected sub-models: ", msg_names, ".")
  }

  for (i in names(object[["submodels"]])) {
    if (is.null(object[["submodels"]][[i]][[1]][["posterior"]])) {
      stop("Sub-model ", i, " has not been estimated yet, so the global model ",
           "cannot be solved.")
    }
  }

  # Last observations
  end_i <- unlist(lapply(object[["submodels"]],
                         function(x){stats::tsp(x[[1]][["data"]][["train"]][["y"]])[2]}))
  end <- unique(end_i)
  if (length(end) > 1) {
    stop("Last date of last period differs across sub-models.\n",
         "Feel free to send a feature request.")
  }

  # Period of the solution ----
  #
  # Defaults to the last period of the estimation sample. Compared with a
  # tolerance rather than for equality, since the times of a 'ts' object are
  # accumulated in floating point and a quarter of a year is not exact.
  times <- stats::time(object[["global"]][["endogen"]])

  if (is.null(period)) {
    period <- which.min(abs(times - end))
    if (abs(times[period] - end) > 1e-6) {
      stop("The last period of the sub-models is not part of the global data.")
    }
  } else {
    if (length(period) != 1 || period < 1 || period > length(times)) {
      stop("Argument 'period' must be a single integer between 1 and ",
           length(times), ".")
    }
  }

  # Sub-models whose coefficients or error covariance matrices move are reduced
  # to the draws of that period, which leaves every sub-model with constant
  # draws and lets the stacking below treat them alike. This comes first, so
  # that only the draws of that period are collected and transformed.
  #
  # Sub-models are stacked in levels, so one in error correction form is
  # rewritten as the level VARX model it implies next. Done here rather than
  # in a step of its own, so that a GVEC model is solved like a GVAR model and
  # a model whose sub-models were chosen from VARX and VECX candidates, and so
  # come in both forms, is solved at all. The global data, the index and the
  # weight matrices describe the units and their weights, which neither form
  # changes.
  # The cointegrating vectors are kept here, because vec_to_var() drops them:
  # they say nothing about the level form and everything about the long-run
  # relations, which the solved model is asked for further down. One entry per
  # sub-model that has any, of the period this model is solved at.
  coint_vectors <- list()

  for (i in names(object[["submodels"]])) {
    model_i <- .submodel_at_period(object[["submodels"]][[i]][[1]],
                                   time = times[period], submodel = i)

    if (inherits(model_i, "vecxsubmodel")) {
      rank_i <- model_i[["model"]][["rank"]]
      beta_i <- model_i[["posterior"]][["beta"]][["coeffs"]]
      if (!is.null(rank_i) && rank_i > 0 && !is.null(beta_i)) {
        coint_vectors[[i]] <- list(rank = rank_i, coeffs = as.matrix(beta_i))
      }
      model_i <- vec_to_var(model_i)
    }
    object[["submodels"]][[i]][[1]] <- model_i
  }
  
  # Endogenous variables
  endogen_i <- lapply(object[["submodels"]], function(x){x[[1]][["model"]][["endogen"]]})
  pos_i <- endogen_i
  for (i in submodels) {
    index_i <- index[index[, "submodel"] == i & index[, "variable"] %in% endogen_i[[i]], ]
    # match() the other way round would give, for every row of index_i, its
    # position among the model's variables -- not the row of index_i that
    # belongs to the j-th variable of the model, which is what is wanted here.
    pos_i[[i]] <- index_i[match(endogen_i[[i]], index_i[, "variable"]), "id"]
    rm(index_i)
    endogen_i[[i]] <- paste0(i, "_", endogen_i[[i]])
  }
  
  # Number of endogenous variables
  k_endogen_i <- unlist(lapply(object[["submodels"]],
                               function(x){x[[1]][["model"]][["k_endogen"]]}))
  # Maximum lag of endogenous variables
  p_endogen_i <- unlist(lapply(object[["submodels"]],
                               function(x){x[[1]][["model"]][["p_endogen"]]}))
  n_endogen_i <- k_endogen_i * k_endogen_i * p_endogen_i
  
  # Number of weakly exogenous variables
  k_exogen_i <- unlist(lapply(object[["submodels"]],
                              function(x){x[[1]][["model"]][["k_exogen"]]}))
  # Maximum lag of exogenous variables
  p_exogen_i <- unlist(lapply(object[["submodels"]],
                              function(x){x[[1]][["model"]][["p_exogen"]]}))
  # 'p_exogen' counts blocks from the contemporaneous term onwards, so a
  # sub-model with p_exogen blocks reaches back to lag p_exogen - 1.
  n_exogen_i <- k_endogen_i * k_exogen_i * p_exogen_i
  
  # Final dimensions of the global model
  k <- sum(k_endogen_i)
  p <- max(p_endogen_i, p_exogen_i - 1, 0)

  # Without a lag nothing propagates between the units, so there is no dynamic
  # global model to solve.
  if (p == 0) {
    stop("None of the sub-models uses a lag of its endogenous or of its weakly ",
         "exogenous variables, so the global model has no dynamics.")
  }
  # Position of the first variable of each sub-model in the global vector.
  # Written without indexing so that a global model of a single sub-model,
  # where the shift would be an empty assignment, works as well.
  k_pos <- cumsum(c(0, k_endogen_i[-length(k_endogen_i)]))
  names(k_pos) <- names(k_endogen_i)
  
  # Use global variables
  m_i <- unlist(lapply(object[["submodels"]],
                       function(x){x[[1]][["model"]][["m_global"]]}))
  m <- unique(m_i)
  if (length(m) > 1) {
    stop("Number of global variables differs across sub-models.\n",
         "Feel free to send a feature request.")
  }
  s_i <- rep(0, length(m_i))
  if (m > 0) {
    s_i <- unlist(lapply(object[["submodels"]],
                         function(x){x[[1]][["model"]][["s_global"]]}))
  }
  # Global variables that contribute no regressor to any sub-model leave no
  # coefficient to carry over, so the global model does not use them either.
  global <- m > 0 & max(s_i) > 0
  if (global) {
    # The highest lag a sub-model reaches is one below its number of blocks.
    s <- max(s_i) - 1
    n_global_i <- k_endogen_i * m_i * s_i
    global_i <- unique(unlist(lapply(object[["submodels"]], function(x){x[[1]][["model"]][["global"]]})))
  } else {
    n_global_i <- k_endogen_i * 0
    s <- 0
    global_i <- NULL
  }
  
  # Deterministic variables
  #
  # The sub-models need not carry the same terms. A VECX sub-model at rank zero
  # has no error correction term, so a trend restricted to it drops out of its
  # level representation, while the sub-models of higher rank keep it. The
  # global model uses every term that any sub-model uses, and a sub-model gets
  # zero coefficients on the terms it lacks.
  n_i <- vapply(object[["submodels"]], function(x) {
    n_x <- x[[1]][["model"]][["n"]]
    if (is.null(n_x)) 0 else as.numeric(n_x)
  }, numeric(1))
  det_names_i <- lapply(object[["submodels"]], function(x){x[[1]][["model"]][["deterministic"]]})
  for (i in submodels) {
    if (length(det_names_i[[i]]) != n_i[i]) {
      stop("Sub-model ", i, " uses ", n_i[i], " deterministic term(s), but ",
           "names ", length(det_names_i[[i]]), " of them.")
    }
  }
  det_names <- unique(unlist(det_names_i))
  names(det_names) <- NULL
  # In the order of the global model's deterministic series, so the result does
  # not depend on which sub-model happens to come first. Terms it does not
  # provide are reported further down, when the data are taken from it.
  det_global <- dimnames(object[["global"]][["deterministic"]])[[2]]
  det_names <- c(intersect(det_global, det_names), setdiff(det_names, det_global))
  n <- length(det_names)
  if (n == 0) {
    det_names <- NULL
  }
  
  # Number of draws
  draws_i <- unlist(lapply(object[["submodels"]],
                           function(x){nrow(x[[1]][["posterior"]][["u_sigma_inv"]][["coeffs"]])}))
  draws <- unique(draws_i)
  if (length(draws) > 1) {
    stop("Number of posterior draws differs across sub-models.\n",
         "Feel free to send a feature request.")
  }
  
  # Labels of the draws. Row i of every sub-model is stacked into draw i of the
  # global model, so the rows have to stand for the same draws of their chains.
  # Sub-models thinned by different intervals, or kept from different stretches
  # of their chains, would otherwise be combined without notice when they happen
  # to keep as many draws.
  mcpar_i <- lapply(object[["submodels"]], function(x) {
    labels <- attr(x[[1]][["posterior"]][["u_sigma_inv"]][["coeffs"]], "mcpar")
    if (is.null(labels)) c(1, draws, 1) else as.numeric(labels)
  })
  mcpar <- mcpar_i[[1]]
  if (!all(vapply(mcpar_i, function(x) isTRUE(all.equal(x, mcpar)), logical(1)))) {
    stop("The posterior draws of the sub-models are labelled differently, so a ",
         "row does not stand for the same draw in all of them: ",
         paste0(names(mcpar_i), " from ", vapply(mcpar_i, `[`, numeric(1), 1),
                " to ", vapply(mcpar_i, `[`, numeric(1), 2), " by ",
                vapply(mcpar_i, `[`, numeric(1), 3), collapse = ", "),
         ". Thin the sub-models alike before solving the global model.")
  }

  # Number of observations
  tt_i <- unlist(lapply(object[["submodels"]],
                        function(x){nrow(x[[1]][["data"]][["train"]][["y"]])}))
  tt <- unique(tt_i)
  if (length(tt) > 1) {
    stop("Number of available observations differs across sub-models.\n",
         "Feel free to send a feature request.")
  }

  # Create the bvarmodel object ----
  
  # Model
  model <- NULL
  model[["type"]] <- "GVAR"
  model[["k"]] <- k
  model[["p"]] <- p
  model[["m"]] <- m
  model[["s"]] <- s
  model[["n"]] <- n
  model[["varsel"]] <- "none"
  endogen <- unlist(endogen_i)
  names(endogen) <- NULL
  model[["endogen"]] <- endogen
  if (global) {
    model[["exogen"]] <- global_i
  }
  model[["deterministic"]] <- det_names
  # The global model is returned in its reduced form (see below), whose errors
  # have a constant covariance matrix.
  model[["error"]] <- "wishart"
  model[["tvp"]] <- FALSE
  # The contemporaneous matrix G is inverted out below, so what is returned is
  # a reduced form model. G itself is a dense matrix, and the structural form
  # of 'bvartools' stores a unit lower triangular matrix of k * (k - 1) / 2
  # free elements, which cannot hold it.
  model[["structural"]] <- FALSE
  model[["iterations"]] <- as.integer(draws)
  model[["burnin"]] <- 0L
  # Carried only when the sub-models were thinned, as bvartools does, so that
  # what the methods of bvartools derive from the draws is labelled alike.
  if (mcpar[3] > 1) {
    model[["thin"]] <- as.integer(mcpar[3])
  }
  # The name every 'bvartools' method that works on draws dispatches on. It
  # describes the three lines above rather than how any sub-model was
  # estimated: the solved global model has constant coefficients, a full error
  # precision and no structural form, whatever the sub-models had, and that
  # combination is what 'bvartools' calls VarNormalWishart. The sub-models are
  # free to differ: a sub-model with time varying parameters or stochastic
  # volatility was reduced to the draws of a single period further up, so it
  # leaves one set of coefficients and one precision matrix per draw as well.
  #
  # Without it add_posterior_forecasts() has nothing to dispatch on and fails
  # with an argument of length zero, which is why a solved global model could
  # not be forecast at all.
  model[["algorithm"]] <- "VarNormalWishart"
  
  # Data
  endogen_data <- object[["global"]][["endogen"]][, unlist(pos_i)]
  y_tsp <- stats::tsp(endogen_data)
  data <- stats::embed(endogen_data, p + 1)
  
  y <- stats::ts(data[, 1:k], end = y_tsp[2], frequency = y_tsp[3])
  dimnames(y)[[2]] <- endogen
  
  x <- stats::ts(data[, -(1:k)], end = y_tsp[2], frequency = y_tsp[3])
  x_names <- NULL
  for (i in 1:p) {
    x_names <- append(x_names, paste0(endogen, ".l", i))
  }
  
  if (global) {
    tsp_global <- stats::tsp(object[["global"]][["exogen"]])
    # The series themselves, before the lags are laid out beside them. This is
    # what 'bvartools' means by data$original$exogen -- one column per variable
    # rather than one per regressor, the way data$original$endogen above is the
    # series and not the lag blocks built from them. prepare_forecast_input()
    # lags whatever it is given, so handing it the embedded matrix would lag it
    # a second time; and it is the only description of the global variables the
    # returned object carries, so a caller assembling their future path has
    # nothing else to read.
    global_original <- stats::ts(
      matrix(object[["global"]][["exogen"]][, global_i], ncol = length(global_i)),
      end = tsp_global[2], frequency = tsp_global[3])
    dimnames(global_original)[[2]] <- global_i

    global_data <- stats::embed(object[["global"]][["exogen"]][, global_i], s + 1)
    global_data <- stats::ts(global_data, end = tsp_global[2], frequency = tsp_global[3])
    x <- stats::na.omit(cbind(x, global_data))
    for (i in 1:(s + 1)) {
      x_names <- append(x_names, paste0(global_i, ".l", i - 1)) 
    }
  } else {
    global_data <- NULL
    global_original <- NULL
  }

  # The global model has to describe the sample the sub-models were estimated
  # on. Embedding the global series only drops the leading observations that
  # the lag structure of the global model needs, which can leave one more than
  # align_model_obs left the sub-models with -- for instance when the models
  # that were selected use fewer lags than the longest candidate did.
  sub_tsp <- stats::tsp(object[["submodels"]][[submodels[1]]][[1]][["data"]][["train"]][["y"]])
  if (sub_tsp[1] < stats::tsp(y)[1] - 1e-6) {
    stop("The sub-models were estimated from ", sub_tsp[1], " on, which the ",
         "lag structure of the global model leaves no room for.")
  }
  y <- stats::window(y, start = sub_tsp[1], end = sub_tsp[2])
  x <- stats::window(x, start = sub_tsp[1], end = sub_tsp[2])

  if (n > 0) {

    # From the global model, which is where the sub-models took them from as
    # well, so the global model uses the same series they were estimated with.
    det_data <- .submodel_deterministic(object, det_names,
                                        start = sub_tsp[1], end = sub_tsp[2])

    if (nrow(det_data) != nrow(x)) {
      stop("The deterministic terms span ", nrow(det_data), " observations, ",
           "the global model ", nrow(x), ".")
    }

    x <- cbind(x, det_data)
    x_names <- append(x_names, det_names)

  } else {
    det_data <- NULL
  }
  dimnames(x)[[2]] <- x_names
  
  data <- list("original" = list("endogen" = endogen_data,
                                 "exogen" = global_original,
                                 "deterministic" = det_data),
               "train" = list("y" = y,
                              "x" = x))
  
  # Create skeleton
  posterior <- NULL
  posterior[["a0"]] <- matrix(NA_real_, draws, k * k)
  posterior[["a"]] <- matrix(NA_real_, draws, k * k * p)
  
  if (global) {
    posterior[["b"]] <- matrix(NA, draws, k * m * (1 + s))
  } else {
    posterior[["b"]] <- NULL
  }
  
  if (n > 0) {
    posterior[["c"]] <- matrix(NA, draws, k * n)
  } else {
    posterior[["c"]] <- NULL
  }
  
  u_sigma_inv <- matrix(NA, draws, k * k)

  # Weight matrices ----
  #
  # Constant weights are stored as a single matrix, which is then used for
  # every period.
  n_periods_w <- vapply(submodels, function(i) {
    n_vars <- length(unique(index[index[, "submodel"] == i, "variable"])) +
      length(unique(index[index[, "submodel"] != i, "variable"]))
    nrow(object[["weights"]][[i]]) / n_vars
  }, numeric(1))
  period_w <- if (all(n_periods_w == 1)) 1 else period

  # Weight matrix ----
  #
  # Rows in the order the coefficients of a sub-model are in, columns
  # restricted to the variables that are endogenous to the global model.
  keep <- unlist(pos_i)
  names(keep) <- NULL
  w <- lapply(submodels, .submodel_weight_matrix, object = object,
              period = period_w, keep = keep)
  names(w) <- submodels
  
  message("Combining sub-models to the global model ...")

  # Only interactively: a progress bar in the output of R CMD check or of a
  # batch job is noise.
  pb <- NULL
  if (interactive()) {
    pb <- utils::txtProgressBar(style = 3)
    on.exit(close(pb), add = TRUE)
  }
  for (draw in 1:draws) {
    
    # Put together A0 ----
    a0_temp <- matrix(NA_real_, k, k)
    for (i in submodels) {
      
      pos_a0_exogen <- n_endogen_i[i] + 1:(k_endogen_i[i] * k_exogen_i[i])
      
      # A0
      if (isTRUE(object[["submodels"]][[i]][[1]][["model"]][["structural"]]) & k_endogen_i[i] > 1) {
        # A structural sub-model carries the k (k - 1) / 2 free elements of its
        # unit lower triangular A0 as the last coefficients. They are stored
        # column by column -- (2,1), (3,1), ..., (k,1), (3,2), ... -- which is
        # the order of which(lower.tri()). The contemporaneous regressors enter
        # z with a negative sign, so the free elements are the entries of A0 in
        # A0 y_t = ... + u_t as they are.
        a_i <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]]
        n_struct <- k_endogen_i[i] * (k_endogen_i[i] - 1) / 2
        A0 <- diag(1, k_endogen_i[i])
        A0[lower.tri(A0)] <- a_i[draw, ncol(a_i) - n_struct + 1:n_struct]
        rm(a_i)
      } else {
        A0 <- diag(1, k_endogen_i[i])
      }
      
      # Contemporary foreign
      A0_exogen <- matrix(object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_a0_exogen], k_endogen_i[i])
      
      a0_temp[k_pos[i] + 1:k_endogen_i[i],] <- cbind(A0, -A0_exogen) %*% w[[i]]
    }
    # A sub-model contributes its own variables through an identity block of
    # its weight matrix, so the diagonal must come out as one. Compared with a
    # tolerance, since it is the result of a matrix product.
    if (any(abs(diag(a0_temp) - 1) > 1e-8)) {
      stop("Not all diagonal elements of the global A0 matrix are 1.")
    }
    posterior[["a0"]][draw, ] <- matrix(a0_temp)
    rm(a0_temp)
    
    # Put together A ----
    for (j in 1:p) {
      
      # Create global matrix of [A_d, A_*] for lag j
      # Fill with NA to facilitate error spotting
      a_temp <- matrix(NA_real_, k, k)
      
      for (i in submodels) {
        
        # Create a sub-model matrix [A_d, A_*] and fill it with A_d and A_*
        temp_i <- matrix(0, k_endogen_i[i], k_endogen_i[i] + k_exogen_i[i])
        
        # Endogenous draws of lag j
        if (j <= p_endogen_i[i]) { # If j is larger than p_endogen_i, leave the country matrix 0
          pos_endogen <- (j - 1) * k_endogen_i[i] * k_endogen_i[i] + 1:(k_endogen_i[i] * k_endogen_i[i])
          temp_i[, 1:k_endogen_i[i]] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_endogen]
        }
        
        # Foreign draws of lag j
        if (j <= p_exogen_i[i] - 1) { # Beyond the blocks of sub-model i, leave the country matrix 0
          pos_exogen <-  n_endogen_i[i] + j * (k_endogen_i[i] * k_exogen_i[i]) + 1:(k_endogen_i[i] * k_exogen_i[i])
          temp_i[, k_endogen_i[i] + 1:k_exogen_i[i]] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_exogen]
        }
        
        a_temp[k_pos[i] + 1:k_endogen_i[i],] <- temp_i %*% w[[i]]
        rm(temp_i)
      }
      
      # Store
      posterior[["a"]][draw, (j - 1) * k^2 + 1:(k^2)] <- matrix(a_temp)
      rm(a_temp)
    }
    
    # Put together H ----
    if (global) {
      for (j in 1:(s + 1)) {
        h_temp <- matrix(0, k, m)
        for (i in submodels) {
          if (m_i[i] > 0) {
            if (j <= s_i[i]) {
              pos_global <- n_endogen_i[i] + n_exogen_i[i] + (j - 1) * k_endogen_i[i] * m_i[i] + 1:(k_endogen_i[i] * m_i[i])
              h_temp[k_pos[i] + 1:k_endogen_i[i],] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_global]
            }  
          }
        }
        
        # Store
        posterior[["b"]][draw, (j - 1) * k * m + 1:(k * m)] <- matrix(h_temp)
        rm(h_temp)
      }
    }
    
    #### Put together D ####
    if (n > 0) {
      # Terms a sub-model does not use keep their coefficients of zero.
      d_temp <- matrix(0, k, n)
      for (i in submodels) {
        if (n_i[i] > 0) {
          pos_det <- n_endogen_i[i] + n_exogen_i[i] + n_global_i[i] + 1:(k_endogen_i[i] * n_i[i])
          # Matched by name, since a sub-model may order its terms differently.
          col_det <- match(det_names_i[[i]], det_names)
          d_temp[k_pos[i] + 1:k_endogen_i[i], col_det] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_det]
        }
      }
      
      # Premultiply by A0_i and store
      posterior[["c"]][draw, ] <- matrix(d_temp)
      rm(d_temp)
    }
    
    # Put together Sigma ----
    sigma_temp <- matrix(0, k, k)
    for (i in submodels) {
      sigma_temp[k_pos[i] + 1:k_endogen_i[i], k_pos[i] + 1:k_endogen_i[i]] <- object[["submodels"]][[i]][[1]][["posterior"]][["u_sigma_inv"]][["coeffs"]][draw, ]
    }
    u_sigma_inv[draw, ] <- sigma_temp
    rm(sigma_temp)
    
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, value = draw / draws)
    }
  }
  
  # Reduced form ----
  #
  # Stacking the sub-models gives the structural form
  #   G y_t = A_1 y_{t-1} + ... + B x_t + C d_t + u_t,
  # in which G is the dense matrix of contemporaneous coefficients assembled
  # above. Every method that follows works on the reduced form
  #   y_t = G^-1 A_1 y_{t-1} + ... + G^-1 B x_t + G^-1 C d_t + G^-1 u_t,
  # so the model is solved for it here rather than leaving the inversion to
  # each of them.
  #
  # G is not triangular, so it cannot be passed on as the structural form of a
  # 'bvarmodel', which stores the k * (k - 1) / 2 free elements of a unit lower
  # triangular matrix. Its draws are returned separately, as element 'g'.
  for (draw in 1:draws) {

    g_i <- matrix(posterior[["a0"]][draw, ], k)
    g_inv <- tryCatch(solve(g_i), error = function(e) {
      stop("The global matrix of contemporaneous coefficients is singular in ",
           "draw ", draw, ", so the global model cannot be solved.")
    })

    for (j in seq_len(p)) {
      pos <- (j - 1) * k^2 + 1:(k^2)
      posterior[["a"]][draw, pos] <- g_inv %*% matrix(posterior[["a"]][draw, pos], k)
    }

    if (global) {
      posterior[["b"]][draw, ] <- g_inv %*% matrix(posterior[["b"]][draw, ], k)
    }

    if (n > 0) {
      posterior[["c"]][draw, ] <- g_inv %*% matrix(posterior[["c"]][draw, ], k)
    }

    # The reduced form error is G^-1 u_t, so its covariance is
    # G^-1 Sigma_u G^-1', and the precision that is stored is G' Sigma_u^-1 G.
    # Written in terms of the precision, so that the block diagonal Sigma_u
    # never has to be inverted.
    omega <- matrix(u_sigma_inv[draw, ], k)
    u_sigma_inv[draw, ] <- t(g_i) %*% omega %*% g_i
  }

  # Posterior. The draws keep the labels of the sub-models they were stacked
  # from.
  as_draws <- function(x) {
    coda::mcmc(x, start = mcpar[1], end = mcpar[2], thin = mcpar[3])
  }
  coeffs <- as_draws(cbind(posterior[["a"]], posterior[["b"]], posterior[["c"]]))
  u_sigma_inv <- as_draws(u_sigma_inv)

  # Cointegrating vectors in global coordinates ----
  #
  # A sub-model's cointegrating vector spans its own variables and its foreign
  # ones, and the weight matrix is what turns that pair into a statement about
  # the global vector -- the same mapping the coefficients above go through,
  # transposed, since beta multiplies the variables rather than being multiplied
  # by them. Stacking them gives the sum(r_i) relations of the global system,
  # which is what a persistence profile follows and what the count of common
  # stochastic trends, k - sum(r_i), is about.
  #
  # Kept out of the draw loop: it needs nothing the stacking produces, only the
  # weight matrices and each sub-model's own draws.
  beta_global <- .global_beta(coint_vectors, submodels, w, k, draws)

  result <- list("model" = model,
                 "data" = data,
                 "posterior" = list("a" = list("coeffs" = coeffs),
                                    "u_sigma_inv" = list("coeffs" = u_sigma_inv)),
                 "g" = as_draws(posterior[["a0"]]))

  class(result) <- list("bvarmodel", "list")
  if (!is.null(beta_global)) {
    labelled <- as_draws(beta_global)
    # coda::mcmc() keeps only what it knows about, so the names of the relations
    # are put back on the object the caller sees. The count travels in
    # model$rank rather than as an attribute: an attribute would survive on a
    # model solved in one piece and be lost on one solved in several, where the
    # pieces are rbound, and the two would then differ by a stray attribute.
    attr(labelled, "relations") <- attr(beta_global, "relations")
    attr(labelled, "rank") <- NULL
    result[["posterior"]][["beta"]] <- list("coeffs" = labelled)
    result[["model"]][["rank"]] <- attr(beta_global, "rank")
    # Prepended, so that every method written for the level VAR keeps working
    # and only what needs the cointegrating vectors dispatches differently.
    #
    # A sub-model already in level form contributes none, because vec_to_var()
    # dropped them and nothing can bring them back. A model whose sub-models
    # come in both forms therefore carries fewer relations than its ranks
    # suggest, and 'rank' counts the ones it actually has -- which is the number
    # the count of common stochastic trends, k - rank, must be taken against.
    class(result) <- list("gvecsolved", "bvarmodel", "list")
  }
  return(result)
}

# The sub-model reduced to its draws at a single date.
#
# A time varying sub-model stores its coefficients period after period, one
# block of ncol(z) per period, and a sub-model with stochastic volatility the
# k x k error precision of every period. A time varying sub-model in error
# correction form stores its cointegration vectors period after period as well.
# What is kept is the block of the period at 'time', so that the result has the
# layout of a sub-model with constant coefficients and a constant error
# covariance matrix. A sub-model that has neither is returned as it is.
.submodel_at_period <- function(model, time, submodel) {

  k <- model[["model"]][["k"]]
  kk <- k * k
  y <- model[["data"]][["train"]][["y"]]
  tt <- nrow(y)
  tvp <- isTRUE(model[["model"]][["tvp"]])
  a <- model[["posterior"]][["a"]][["coeffs"]]
  u_sigma_inv <- model[["posterior"]][["u_sigma_inv"]][["coeffs"]]
  sv <- ncol(u_sigma_inv) == tt * kk && tt > 1

  if (!tvp && !sv) {
    return(model)
  }

  # The block of a period keeps the labels of the draws it is cut from, so that
  # a thinned chain still says which draws it kept.
  keep_labels <- function(block, draws) {
    labels <- attr(draws, "mcpar")
    if (is.null(labels)) {
      return(coda::mcmc(block))
    }
    coda::mcmc(block, start = labels[1], end = labels[2], thin = labels[3])
  }

  tsp_y <- stats::tsp(y)
  pos <- round((time - tsp_y[1]) * tsp_y[3]) + 1
  if (pos < 1 || pos > tt) {
    stop("Sub-model ", submodel, " has time varying parameters or stochastic ",
         "volatility, which are only available for its estimation sample from ",
         tsp_y[1], " to ", tsp_y[2], ", so the model cannot be solved at ",
         time, ".")
  }

  if (tvp) {
    n_a <- ncol(a) / tt
    if (n_a != round(n_a)) {
      stop("The coefficient draws of sub-model ", submodel, " do not divide ",
           "into ", tt, " periods.")
    }
    model[["posterior"]][["a"]] <- list("coeffs" = keep_labels(a[, (pos - 1) * n_a + 1:n_a, drop = FALSE], a))
    beta <- model[["posterior"]][["beta"]][["coeffs"]]
    if (!is.null(beta)) {
      n_beta <- ncol(beta) / tt
      model[["posterior"]][["beta"]] <- list("coeffs" = keep_labels(beta[, (pos - 1) * n_beta + 1:n_beta, drop = FALSE], beta))
    }
    model[["posterior"]][["psi"]] <- NULL
    model[["model"]][["tvp"]] <- FALSE
  }

  if (sv) {
    model[["posterior"]][["u_sigma_inv"]] <- list("coeffs" = keep_labels(u_sigma_inv[, (pos - 1) * kk + 1:kk, drop = FALSE], u_sigma_inv))
    model[["posterior"]][["u_omega_inv"]] <- NULL
    model[["model"]][["error"]] <- sub("sv", "gamma", model[["model"]][["error"]], fixed = TRUE)
  }

  return(model)
}

# The cointegrating vectors of the sub-models, in the coordinates of the global
# model.
#
# 'coint_vectors' is what .combine_submodels() kept before vec_to_var() dropped
# it, one entry per sub-model that has a long-run relation.
#
# 'w[[i]]' has one row per variable of sub-model i, its own first and then its
# foreign ones, and one column per variable of the global model. A vector of the
# sub-model is stated over the same variables, in the same order, followed by
# any restricted deterministic term, which is dropped: a deterministic term is
# not a variable of the global system and the profile of a relation is about the
# variables in it.
#
# Returns draws x (k * sum(r_i)), the relations side by side, or NULL where no
# sub-model has one -- every sub-model at rank zero, or a global VAR.
.global_beta <- function(coint_vectors, submodels, w, k, draws) {

  blocks <- list()
  labels <- character(0)

  for (i in submodels) {
    kept <- coint_vectors[[i]]
    if (is.null(kept)) {
      next
    }
    rank_i <- kept[["rank"]]
    coeffs <- kept[["coeffs"]]

    mapped <- nrow(w[[i]])
    n_beta <- ncol(coeffs) / rank_i
    if (n_beta < mapped) {
      stop("Sub-model ", i, " has ", n_beta, " rows in its cointegrating ",
           "vectors but its weight matrix maps ", mapped, " variables.")
    }

    block <- matrix(NA_real_, draws, k * rank_i)
    for (draw in seq_len(draws)) {
      beta_i <- matrix(coeffs[draw, ], n_beta, rank_i)
      block[draw, ] <- as.numeric(crossprod(w[[i]], beta_i[seq_len(mapped), , drop = FALSE]))
    }
    blocks[[length(blocks) + 1]] <- block
    labels <- c(labels, paste0(i, "_r", seq_len(rank_i)))
  }

  if (length(blocks) == 0) {
    return(NULL)
  }

  out <- do.call(cbind, blocks)
  attr(out, "rank") <- length(labels)
  attr(out, "relations") <- labels
  out
}
