test_that("missing or malformed elements are rejected", {
  submodel_data <- gvar_data()

  no_endogen <- submodel_data
  no_endogen[["US"]][["endogen"]] <- NULL
  expect_error(bgvars:::.check_submodeldata(no_endogen),
               "does not contain element 'endogen'")

  no_ts_endogen <- submodel_data
  no_ts_endogen[["US"]][["endogen"]] <- unclass(submodel_data[["US"]][["endogen"]])
  expect_error(bgvars:::.check_submodeldata(no_ts_endogen),
               "must be a time-series object")

  no_weights <- submodel_data
  no_weights[["JP"]][["weights"]] <- NULL
  expect_error(bgvars:::.check_submodeldata(no_weights),
               "does not contain element 'weights'")

  no_ts_weights <- submodel_data
  no_ts_weights[["JP"]][["weights"]] <- unclass(submodel_data[["JP"]][["weights"]])
  expect_error(bgvars:::.check_submodeldata(no_ts_weights),
               "must be a time-series object")
})

test_that("non-zero own weights are rejected", {
  submodel_data <- gvar_data()
  submodel_data[["CA"]][["weights"]][, "CA"] <- 1

  expect_error(bgvars:::.check_submodeldata(submodel_data),
               "own weights must be zero")
})
