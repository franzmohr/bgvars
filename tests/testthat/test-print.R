test_that("print.submodeldata shows the variables of every sub-model", {
  submodel_data <- gvar_data()

  output <- capture.output(print(submodel_data))
  expect_match(output[1], "submodel")
  expect_length(output, length(submodel_data) + 1L)
  for (i in names(submodel_data)) {
    expect_true(any(grepl(i, output, fixed = TRUE)))
  }
  # Every variable of the data set becomes a column.
  for (i in dimnames(submodel_data[["US"]][["endogen"]])[[2]]) {
    expect_true(any(grepl(i, output, fixed = TRUE)))
  }
})

test_that("print.submodeldata returns the availability table invisibly", {
  submodel_data <- gvar_data()

  utils::capture.output(result <- expect_invisible(print(submodel_data)))
  expect_s3_class(result, "data.frame")
  expect_equal(result[, "submodel"], names(submodel_data))
})

test_that("print.submodeldata marks variables that are missing for a sub-model", {
  submodel_data <- gvar_data()
  submodel_data[["JP"]][["endogen"]] <-
    submodel_data[["JP"]][["endogen"]][, c("y", "Dp")]

  utils::capture.output(table <- print(submodel_data))
  expect_equal(table[table[, "submodel"] == "JP", "eq"], "")
  expect_equal(table[table[, "submodel"] == "US", "eq"], "x")
})

test_that("print.gvarmodel reports variables and global data", {
  object <- gvar_object()

  output <- capture.output(print(object))
  expect_true(any(grepl("Vector Autoregressive Model", output)))
  expect_true(any(grepl("Variables in sub-models", output)))
  expect_true(any(grepl("Global variables", output)))
  expect_true(any(grepl("poil", output, fixed = TRUE)))
})

test_that("print.gvarmodel omits the global block if there are no global data", {
  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data)

  output <- capture.output(print(object))
  expect_false(any(grepl("Global variables", output)))
})
