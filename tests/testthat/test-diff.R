test_that("diff.submodeldata differences the selected variables only", {
  submodel_data <- gvar_data()
  result <- diff(submodel_data, variables = "y")

  original <- submodel_data[["US"]][["endogen"]]
  differenced <- result[["US"]][["endogen"]]

  expect_equal(names(result), names(submodel_data))
  expect_equal(nrow(differenced), nrow(original) - 1L)
  expect_equal(dimnames(differenced)[[2]], dimnames(original)[[2]])
  expect_equal(as.numeric(differenced[, "y"]), as.numeric(diff(original[, "y"])))
  # Variables that were not selected are only shortened, not differenced.
  expect_equal(as.numeric(differenced[, "Dp"]), as.numeric(original[-1, "Dp"]))
})

test_that("diff.submodeldata applies the multiplier", {
  submodel_data <- gvar_data()
  result <- diff(submodel_data, variables = "y", multi = 100)

  original <- submodel_data[["US"]][["endogen"]]
  expect_equal(as.numeric(result[["US"]][["endogen"]][, "y"]),
               as.numeric(diff(original[, "y"])) * 100)
})

test_that("diff.submodeldata differences all variables by default", {
  submodel_data <- gvar_data()
  result <- diff(submodel_data)

  original <- submodel_data[["JP"]][["endogen"]]
  for (i in dimnames(original)[[2]]) {
    expect_equal(as.numeric(result[["JP"]][["endogen"]][, i]),
                 as.numeric(diff(original[, i])))
  }
})

test_that("diff.submodeldata shifts the start of the series", {
  submodel_data <- gvar_data()
  result <- diff(submodel_data, variables = "y")

  original <- stats::tsp(submodel_data[["US"]][["endogen"]])
  expect_equal(stats::tsp(result[["US"]][["endogen"]]),
               c(original[1] + 1 / original[3], original[2], original[3]))
})

test_that("diff.submodeldata leaves the weight series untouched", {
  submodel_data <- gvar_data()
  result <- diff(submodel_data, variables = "y")

  expect_equal(result[["US"]][["weights"]], submodel_data[["US"]][["weights"]])
})

test_that("diff.submodeldata rejects unknown variables", {
  submodel_data <- gvar_data()
  expect_error(diff(submodel_data, variables = "does_not_exist"),
               "is contained in the data")
})
