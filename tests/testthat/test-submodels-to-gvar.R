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
                                exogen = c("y", "Dp"), p_exogen = 0,
                                global = "poil", s = 0,
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
                               exogen = c("y", "r"), p_exogen = 0,
                               global = "poil", s = 0,
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
                                exogen = c("y", "eq"), p_exogen = 0,
                                global = "poil", s = 0,
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
                                exogen = c("y", "Dp"), p_exogen = 0,
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
                          exogen = c("y", "Dp"), p_exogen = 0,
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

test_that("submodels_to_gvar refuses sub-models in error correction form", {
  object <- gvec_estimated(iterations = 20, burnin = 10)

  expect_error(submodels_to_gvar(object),
               "error correction form")
})

test_that("a global model without any lag is reported", {
  object <- gvar_estimated_spec(endogen = c("y", "Dp"), p_endogen = 0,
                                exogen = c("y", "Dp"), p_exogen = 0,
                                deterministic = "const")

  expect_error(submodels_to_gvar(object), "no dynamics")
})
