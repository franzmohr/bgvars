# Solving the global model.
#
# The tests below pin down two things that used to go wrong silently: the
# variables of the global model have to be the series they are named after, and
# what comes back has to be the reduced form of the stacked sub-models rather
# than something of merely the right shape.

test_that("submodels_to_gvar solves the global model", {
  object <- gvar_estimated(iterations = 30, burnin = 10)

  gvar <- submodels_to_gvar(object)

  expect_s3_class(gvar, "bvarmodel")

  k <- sum(vapply(object[["submodels"]],
                  function(x) x[[1]][["model"]][["k_endogen"]], numeric(1)))
  expect_equal(gvar[["model"]][["k"]], k)
  expect_equal(gvar[["model"]][["type"]], "GVAR")
  expect_equal(length(gvar[["model"]][["endogen"]]), k)
  expect_equal(nrow(gvar[["posterior"]][["a"]][["coeffs"]]), 30L)
})

test_that("the variables of the global model are named after their sub-model", {
  object <- gvar_estimated(iterations = 20, burnin = 10)

  gvar <- submodels_to_gvar(object)

  expected <- unlist(lapply(names(object[["submodels"]]), function(i) {
    paste0(i, "_", object[["submodels"]][[i]][[1]][["model"]][["endogen"]])
  }))
  expect_equal(gvar[["model"]][["endogen"]], expected)
  expect_equal(dimnames(gvar[["data"]][["train"]][["y"]])[[2]], expected)
})

test_that("a variable of the global model carries the series it is named after", {
  # Regression test. The sub-models order their variables the way they were
  # named in add_submodels(), the global index orders them the way the data do.
  # The two were matched the wrong way round, which labelled the columns of the
  # global model with the wrong series whenever the orders disagreed.
  object <- gvar_estimated_spec(endogen = c("Dp", "y"), p_endogen = 1,
                                exogen = c("y", "Dp"), p_exogen = 1,
                                global = "poil", s = 1,
                                deterministic = "const")

  gvar <- submodels_to_gvar(object)

  expect_equal(gvar[["model"]][["endogen"]],
               c("US_Dp", "US_y", "JP_Dp", "JP_y", "CA_Dp", "CA_y"))

  global <- object[["global"]][["endogen"]]
  y <- gvar[["data"]][["train"]][["y"]]
  for (variable in colnames(y)) {
    expect_equal(as.numeric(utils::tail(y[, variable], 5)),
                 as.numeric(utils::tail(global[, variable], 5)),
                 info = variable)
  }
})

test_that("sub-models may use a subset of the available variables", {
  # Regression test. The dimension of the global model follows the endogenous
  # variables of the sub-models, while the weight matrices span every variable
  # in the global index. Taking the weight matrices unrestricted made the two
  # disagree as soon as a sub-model left a variable out.
  object <- gvar_estimated_spec(endogen = c("r", "y"), p_endogen = 1,
                               exogen = c("y", "r"), p_exogen = 1,
                               global = "poil", s = 1,
                               deterministic = "const")

  gvar <- submodels_to_gvar(object)

  expect_equal(gvar[["model"]][["k"]], 6L)
  expect_equal(gvar[["model"]][["endogen"]],
               c("US_r", "US_y", "JP_r", "JP_y", "CA_r", "CA_y"))
  expect_false(any(is.na(gvar[["posterior"]][["a"]][["coeffs"]])))
})

test_that("a weakly exogenous variable no sub-model uses is reported", {
  # "eq" is weakly exogenous to every sub-model but endogenous to none, so the
  # global model cannot be closed and the weight matrices cannot simply be
  # trimmed to the variables that remain.
  object <- gvar_estimated_spec(endogen = c("y", "Dp"), p_endogen = 1,
                                exogen = c("y", "eq"), p_exogen = 1,
                                global = "poil", s = 1,
                                deterministic = "const")

  expect_error(submodels_to_gvar(object),
               "not endogenous to any sub-model")
})

test_that("submodels_to_gvar returns the reduced form", {
  object <- gvar_estimated(iterations = 30, burnin = 10)

  gvar <- submodels_to_gvar(object)

  # A structural 'bvarmodel' stores the free elements of a unit lower
  # triangular matrix, which the dense G of a global model does not fit into.
  expect_false(gvar[["model"]][["structural"]])

  model <- gvar[["model"]]
  k <- model[["k"]]
  expect_equal(ncol(gvar[["posterior"]][["a"]][["coeffs"]]),
               k * k * model[["p"]] +
                 k * model[["m"]] * (model[["s"]] + 1) +
                 k * model[["n"]])

  # The draws of G are kept, and its diagonal is one by construction.
  expect_equal(ncol(gvar[["g"]]), k * k)
  expect_equal(nrow(gvar[["g"]]), nrow(gvar[["posterior"]][["a"]][["coeffs"]]))
  expect_equal(diag(matrix(gvar[["g"]][1, ], k)), rep(1, k))
})

test_that("the reduced form is the structural form premultiplied by G inverse", {
  object <- gvar_estimated(iterations = 20, burnin = 10)

  gvar <- submodels_to_gvar(object)
  k <- gvar[["model"]][["k"]]

  # Rebuild the first lag of the structural form from the sub-models and check
  # that the reduced form is what solving it gives.
  index <- object[["global"]][["index"]]
  submodels <- names(object[["submodels"]])
  keep <- unlist(lapply(submodels, function(i) {
    rows <- index[index[, "submodel"] == i, ]
    rows[match(object[["submodels"]][[i]][[1]][["model"]][["endogen"]],
               rows[, "variable"]), "id"]
  }))

  draw <- 1L
  offset <- 0L
  structural_a <- matrix(NA_real_, k, k)
  for (i in submodels) {
    model <- object[["submodels"]][[i]][[1]][["model"]]
    k_i <- model[["k_endogen"]]
    coeffs <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]]
    a_i <- matrix(coeffs[draw, 1:(k_i * k_i)], k_i)
    w_i <- bgvars:::.submodel_weight_matrix(object, i, period = 1, keep = keep)
    temp <- matrix(0, k_i, k_i + model[["k_exogen"]])
    temp[, 1:k_i] <- a_i
    structural_a[offset + 1:k_i, ] <- temp %*% w_i
    offset <- offset + k_i
  }

  g <- matrix(gvar[["g"]][draw, ], k)
  expect_equal(matrix(gvar[["posterior"]][["a"]][["coeffs"]][draw, 1:(k * k)], k),
               solve(g) %*% structural_a)
})

test_that("generalised impulse responses of the global model are the textbook ones", {
  object <- gvar_estimated(iterations = 30, burnin = 10)
  gvar <- submodels_to_gvar(object)

  n_ahead <- 5L
  result <- irf(gvar, impulse = "US_y", response = "JP_y",
                n_ahead = n_ahead, ci = 0.68, type = "gir", shock = "sd")

  reference <- reference_girf(gvar, impulse = "US_y", response = "JP_y",
                              n_ahead = n_ahead)
  reference <- t(apply(reference, 1, stats::quantile,
                       probs = c(0.16, 0.5, 0.84)))

  expect_equal(as.matrix(result), reference, ignore_attr = TRUE)
})

test_that("a global model without global variables can be solved", {
  object <- gvar_estimated_spec(endogen = c("y", "Dp"), p_endogen = 1,
                                exogen = c("y", "Dp"), p_exogen = 1,
                                deterministic = "const")

  gvar <- submodels_to_gvar(object)

  expect_equal(gvar[["model"]][["m"]], 0L)
  expect_null(gvar[["model"]][["exogen"]])
  expect_equal(ncol(gvar[["posterior"]][["a"]][["coeffs"]]),
               gvar[["model"]][["k"]]^2 * gvar[["model"]][["p"]] +
                 gvar[["model"]][["k"]] * gvar[["model"]][["n"]])
})

test_that("submodels_to_gvar requires one model per sub-model", {
  object <- gvar_estimated_grid(iterations = 20, burnin = 10)

  expect_error(submodels_to_gvar(object),
               "may only contain one model per sub-model")
})

test_that("submodels_to_gvar validates the period of the weight matrices", {
  object <- gvar_estimated(iterations = 20, burnin = 10)

  expect_error(submodels_to_gvar(object, period = 0), "between 1 and")
  expect_error(submodels_to_gvar(object, period = 1e6), "between 1 and")
  expect_error(submodels_to_gvar(object, period = c(1, 2)), "between 1 and")
})

test_that("constant weights make the period irrelevant", {
  # gvar_object() uses a rolling window, so the weights do vary. With weights
  # that were calculated over a fixed set of periods there is one matrix, and
  # every period has to resolve to it.
  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())
  object <- add_weight_matrices(object = object, submodel_data = submodel_data,
                                period = 2014:2016)
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          deterministic = "const",
                          iterations = 20, burnin = 10)
  object <- align_model_obs(object)
  object <- add_priors(object, coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  expect_equal(submodels_to_gvar(object, period = 1)[["g"]],
               submodels_to_gvar(object, period = 20)[["g"]])
})

test_that("candidates in error correction form are refused before they are rewritten", {
  # Only a single model per sub-model is rewritten in levels and stacked, so a
  # GVEC model that still holds several candidates is reported as such rather
  # than converted first.
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", r = 0:1,
                          iterations = 20, burnin = 10)

  expect_error(submodels_to_gvar(object), "only contain one model per sub-model")
})

test_that("a global model without any lag is reported", {
  object <- gvar_estimated_spec(endogen = c("y", "Dp"), p_endogen = 0,
                                exogen = c("y", "Dp"), p_exogen = 1,
                                deterministic = "const")

  expect_error(submodels_to_gvar(object), "no dynamics")
})

# Forecasting the solved global model.
#
# The global model is returned as a 'bvarmodel' so that the methods of
# bvartools work on it, and forecasting is the one that needs more than the
# draws: it dispatches on model$algorithm, which the solved model did not carry
# at all. Without it add_posterior_forecasts() failed on an argument of length
# zero -- the sort of failure that reads as a broken object rather than as a
# missing field.

# The future path of the global variables a forecast needs. Held at the last
# observed value, which is enough to have one and says nothing about the model.
gvar_future_exogen <- function(gvar, n_ahead) {
  exogen <- gvar[["data"]][["original"]][["exogen"]]
  tsp_exogen <- stats::tsp(exogen)
  future <- stats::ts(matrix(rep(as.numeric(exogen[nrow(exogen), ]), each = n_ahead + 1),
                             n_ahead + 1, ncol(exogen)),
                      start = tsp_exogen[2], frequency = tsp_exogen[3])
  dimnames(future)[[2]] <- dimnames(exogen)[[2]]
  future
}

test_that("the solved global model names the algorithm its draws came from", {
  gvar <- submodels_to_gvar(gvar_estimated(iterations = 20, burnin = 10))

  # The reduced form of the stacked sub-models: constant coefficients, a full
  # error precision, no structural form. Whatever the sub-models were.
  expect_equal(gvar[["model"]][["algorithm"]], "VarNormalWishart")
  expect_false(gvar[["model"]][["tvp"]])
  expect_false(gvar[["model"]][["structural"]])
  expect_equal(gvar[["model"]][["error"]], "wishart")
})

test_that("data$original$exogen holds the global series, not their lags", {
  # What prepare_forecast_input() lags to build the forecast regressors, so the
  # embedded matrix that goes into data$train$x would be lagged twice over.
  object <- gvar_estimated(iterations = 20, burnin = 10)
  gvar <- submodels_to_gvar(object)

  exogen <- gvar[["data"]][["original"]][["exogen"]]
  expect_equal(ncol(exogen), gvar[["model"]][["m"]])
  expect_equal(dimnames(exogen)[[2]], gvar[["model"]][["exogen"]])
})

test_that("the solved global model forecasts", {
  gvar <- submodels_to_gvar(gvar_estimated(iterations = 20, burnin = 10))
  k <- gvar[["model"]][["k"]]
  n_ahead <- 4

  gvar <- add_forecast_input(gvar, n_ahead = n_ahead,
                             exogen = gvar_future_exogen(gvar, n_ahead))

  # The compact layout: one row per period, one column per regressor.
  expect_equal(nrow(gvar[["data"]][["forecast"]][["x"]]), n_ahead)
  expect_equal(ncol(gvar[["data"]][["forecast"]][["x"]]),
               k * gvar[["model"]][["p"]] +
                 gvar[["model"]][["m"]] * (gvar[["model"]][["s"]] + 1) +
                 gvar[["model"]][["n"]])

  gvar <- add_posterior_forecasts(gvar)
  forecast <- stats::predict(gvar, n_ahead = n_ahead)

  expect_s3_class(forecast, "bvarprd")
  expect_equal(dim(forecast[["fcst"]]), c(n_ahead, k, 20L))
  expect_true(all(is.finite(forecast[["fcst"]])))
  expect_equal(dimnames(forecast[["fcst"]])[[2]], gvar[["model"]][["endogen"]])
})

test_that("the coefficients of the global model line up with its regressors", {
  # A forecast built on a misread coefficient vector still runs and is still
  # finite -- it is simply a different model. The in-sample residual is what
  # separates the two: against the right layout it is the model's error, and
  # against the same numbers read the other way it is the size of the series.
  gvar <- submodels_to_gvar(gvar_estimated(iterations = 20, burnin = 10))

  k <- gvar[["model"]][["k"]]
  y <- gvar[["data"]][["train"]][["y"]]
  x <- gvar[["data"]][["train"]][["x"]]
  a <- as.matrix(gvar[["posterior"]][["a"]][["coeffs"]])

  expect_equal(ncol(a), k * ncol(x))

  scale <- sqrt(mean(as.matrix(y)^2))
  residual <- function(coefficients) sqrt(mean((t(y) - coefficients %*% t(x))^2))

  # matrix(a[1, ], nrow = k) is the reshape the forecast does.
  expect_lt(residual(matrix(a[1, ], nrow = k)), scale)
  expect_gt(residual(t(matrix(a[1, ], ncol = k))), scale)
})
