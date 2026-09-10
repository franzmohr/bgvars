test_that("create_gvarmodel returns the documented structure", {
  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())

  expect_s3_class(object, "gvarmodel")
  expect_equal(names(object), c("global", "weights", "submodels"))
  expect_equal(names(object[["global"]]),
               c("endogen", "exogen", "deterministic", "index"))

  # The deterministic terms are built once, on the time axis of the global
  # data, so that every sub-model uses the same series.
  deterministic <- object[["global"]][["deterministic"]]
  expect_equal(stats::tsp(deterministic),
               stats::tsp(object[["global"]][["endogen"]]))
  expect_equal(dimnames(deterministic)[[2]],
               c("const", "trend", "season.1", "season.2", "season.3"))
  expect_true(all(deterministic[, "const"] == 1))
  expect_equal(as.numeric(deterministic[, "trend"]), seq_len(nrow(deterministic)))
  expect_null(object[["weights"]])
  expect_null(object[["submodels"]])
})

test_that("the global index maps every sub-model variable exactly once", {
  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())

  index <- object[["global"]][["index"]]
  expect_s3_class(index, "data.frame")
  expect_equal(names(index), c("submodel", "variable", "index", "id"))

  n_vars <- sum(vapply(submodel_data,
                       function(x) NCOL(x[["endogen"]]), numeric(1)))
  expect_equal(nrow(index), n_vars)
  expect_equal(index[, "id"], seq_len(n_vars))
  expect_equal(index[, "index"],
               paste0(index[, "submodel"], "_", index[, "variable"]))
  expect_equal(unique(index[, "submodel"]), names(submodel_data))
})

test_that("global endogenous data are stacked in index order", {
  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data,
                             global_data = gvar_global_data())

  endogen <- object[["global"]][["endogen"]]
  index <- object[["global"]][["index"]]

  expect_s3_class(endogen, "ts")
  expect_equal(ncol(endogen), nrow(index))
  expect_equal(dimnames(endogen)[[2]], index[, "index"])
  expect_equal(stats::tsp(endogen),
               stats::tsp(submodel_data[["US"]][["endogen"]]))
  expect_equal(as.numeric(endogen[, "JP_y"]),
               as.numeric(submodel_data[["JP"]][["endogen"]][, "y"]))
})

test_that("global data are optional and univariate series get a name", {
  submodel_data <- gvar_data()

  without_global <- create_gvarmodel(submodel_data = submodel_data)
  expect_null(without_global[["global"]][["exogen"]])

  poil <- gvar_global_data()[, "poil"]
  expect_null(dimnames(poil))
  univariate <- create_gvarmodel(submodel_data = submodel_data,
                                 global_data = poil)
  expect_equal(dimnames(univariate[["global"]][["exogen"]])[[2]], "global")
  expect_equal(stats::tsp(univariate[["global"]][["exogen"]]), stats::tsp(poil))
})

test_that("create_gvarmodel validates its arguments", {
  submodel_data <- gvar_data()

  expect_error(create_gvarmodel(submodel_data = unclass(submodel_data)),
               "must be of class 'submodeldata'")
  expect_error(create_gvarmodel(submodel_data = submodel_data,
                                global_data = as.data.frame(gvar_global_data())),
               "must be of class 'ts'")
})
