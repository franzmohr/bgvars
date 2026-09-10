test_that("sub-model weights are normalised in every period", {
  submodel_data <- gvar_data()
  weights <- bgvars:::.create_weights_submodel(submodel_data[["US"]], period = 3)

  expect_s3_class(weights, "ts")
  expect_equal(dim(weights),
               c(NROW(submodel_data[["US"]][["endogen"]]),
                 NCOL(submodel_data[["US"]][["weights"]])))
  expect_equal(dimnames(weights)[[2]],
               dimnames(submodel_data[["US"]][["weights"]])[[2]])
  expect_equal(stats::tsp(weights),
               stats::tsp(submodel_data[["US"]][["endogen"]]))
  expect_equal(as.numeric(rowSums(weights)), rep(1, nrow(weights)))
  expect_true(all(weights[, "US"] == 0))
})

test_that("a period vector produces constant weights", {
  submodel_data <- gvec_data()
  weights <- bgvars:::.create_weights_submodel(submodel_data[["US"]],
                                               period = 1999:2001)

  expect_equal(as.numeric(rowSums(weights)), rep(1, nrow(weights)))
  # All rows are identical when weights are not time varying.
  expect_equal(nrow(unique(as.matrix(weights))), 1L)
})

test_that("rolling window weights do vary over time", {
  submodel_data <- gvar_data()
  weights <- bgvars:::.create_weights_submodel(submodel_data[["US"]], period = 3)

  expect_gt(nrow(unique(as.matrix(weights))), 1L)
})

test_that(".create_weights_submodel requires weights of class 'ts'", {
  submodel_data <- gvar_data()
  submodel_data[["US"]][["weights"]] <-
    unclass(submodel_data[["US"]][["weights"]])

  expect_error(bgvars:::.create_weights_submodel(submodel_data[["US"]], period = 3),
               "must be of class 'ts'")
})

test_that("add_weight_matrices adds one matrix per sub-model", {
  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())
  object <- add_weight_matrices(object = object,
                                submodel_data = submodel_data, period = 3)

  expect_equal(names(object[["weights"]]), names(submodel_data))

  index <- object[["global"]][["index"]]
  tt <- nrow(object[["global"]][["endogen"]])
  for (i in names(submodel_data)) {
    n_endogen <- sum(index[, "submodel"] == i)
    n_exogen <- length(unique(index[index[, "submodel"] != i, "variable"]))
    expect_equal(dim(object[["weights"]][[i]]),
                 c((n_endogen + n_exogen) * tt, nrow(index)))
  }
})

test_that("add_weight_matrices validates the sub-model data", {
  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())

  broken <- submodel_data
  broken[["US"]][["weights"]] <- NULL
  expect_error(add_weight_matrices(object, broken, period = 3),
               "does not contain element 'weights'")
})

test_that("get_weight_matrix links a sub-model to the global variable vector", {
  object <- gvar_object()
  w <- get_weight_matrix(object, "US")

  index <- object[["global"]][["index"]]
  vars_endogen <- index[index[, "submodel"] == "US", "variable"]
  vars_exogen <- unique(index[index[, "submodel"] != "US", "variable"])

  expect_equal(dimnames(w)[[1]], c(vars_endogen, vars_exogen))
  expect_equal(dimnames(w)[[2]], index[, "index"])

  # Own variables are picked out by an identity block.
  own <- w[seq_along(vars_endogen), paste0("US_", vars_endogen)]
  expect_equal(own, diag(1, length(vars_endogen)), ignore_attr = TRUE)

  # A sub-model never contributes to its own weakly exogenous variables ...
  foreign <- w[length(vars_endogen) + seq_along(vars_exogen), , drop = FALSE]
  expect_true(all(foreign[, index[, "submodel"] == "US"] == 0))

  # ... and the remaining weights of each weakly exogenous variable add up to one.
  expect_equal(as.numeric(rowSums(foreign)), rep(1, length(vars_exogen)))
})

test_that("weakly exogenous variables are weighted averages of foreign series", {
  object <- gvar_object()
  w <- get_weight_matrix(object, "US")
  endogen <- object[["global"]][["endogen"]]
  index <- object[["global"]][["index"]]

  n_endogen <- sum(index[, "submodel"] == "US")
  vars_exogen <- unique(index[index[, "submodel"] != "US", "variable"])

  # Row of the weakly exogenous variable "y" of the US model. Row names are not
  # unique, since endogenous and weakly exogenous variables share their names.
  row_y <- w[n_endogen + which(vars_exogen == "y"), ]

  expect_equal(as.numeric(row_y %*% endogen[1, ]),
               as.numeric(row_y["JP_y"] * endogen[1, "JP_y"] +
                            row_y["CA_y"] * endogen[1, "CA_y"]))
})

test_that("get_weight_matrix validates its arguments", {
  object <- gvar_object()

  expect_error(get_weight_matrix(object, c("US", "JP")),
               "may only contain one element")
  expect_error(get_weight_matrix(object, "XX"),
               "No weight matrix available")
  expect_error(get_weight_matrix(object, "US", period = 0),
               "must be a single integer")
  expect_error(get_weight_matrix(object, "US", period = 1e6),
               "must be a single integer")

  plain <- create_gvarmodel(submodel_data = gvar_data(),
                            global_data = gvar_global_data())
  expect_error(get_weight_matrix(plain, "US"),
               "does not contain weight matrices")
})

test_that("get_weight_matrix cuts out the matrix of the requested period", {
  object <- gvar_object()

  index <- object[["global"]][["index"]]
  n_vars <- sum(index[, "submodel"] == "US") +
    length(unique(index[index[, "submodel"] != "US", "variable"]))

  for (period in c(1, 5)) {
    expect_equal(get_weight_matrix(object, "US", period = period),
                 object[["weights"]][["US"]][(period - 1) * n_vars + 1:n_vars, ],
                 ignore_attr = TRUE, info = period)
  }

  # The weights are rolling window sums, so different periods really differ.
  expect_false(isTRUE(all.equal(get_weight_matrix(object, "US", period = 1),
                                get_weight_matrix(object, "US", period = 30))))
})

test_that("add_weight_matrices refuses a global model of one sub-model", {
  submodel_data <- gvar_data("US")
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())

  expect_error(add_weight_matrices(object, submodel_data, period = 3),
               "requires at least two sub-models")
})
