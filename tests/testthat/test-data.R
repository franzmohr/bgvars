test_that("shipped data sets have the documented structure", {
  for (name in c("gvar2019", "gvar2023", "dees2007")) {
    utils::data(list = name, package = "bgvars", envir = environment())
    dataset <- get(name, envir = environment())

    expect_type(dataset, "list")
    expect_true("submodel_data" %in% names(dataset))
    expect_true("global_data" %in% names(dataset))
    expect_s3_class(dataset[["submodel_data"]], "submodeldata")
    expect_s3_class(dataset[["global_data"]], "ts")
    expect_false(is.null(dimnames(dataset[["global_data"]])[[2]]))
  }
})

test_that("only the GVAR data sets ship regional weights", {
  utils::data("gvar2019", package = "bgvars", envir = environment())
  utils::data("gvar2023", package = "bgvars", envir = environment())
  utils::data("dees2007", package = "bgvars", envir = environment())

  expect_s3_class(gvar2019[["region_weights"]], "ts")
  expect_s3_class(gvar2023[["region_weights"]], "ts")
  expect_null(dees2007[["region_weights"]])
})

test_that("every sub-model provides named endogenous and weight series", {
  for (name in c("gvar2019", "gvar2023", "dees2007")) {
    utils::data(list = name, package = "bgvars", envir = environment())
    submodel_data <- get(name, envir = environment())[["submodel_data"]]

    expect_false(is.null(names(submodel_data)))
    expect_false(any(duplicated(names(submodel_data))))

    for (i in names(submodel_data)) {
      expect_s3_class(submodel_data[[i]][["endogen"]], "ts")
      expect_s3_class(submodel_data[[i]][["weights"]], "ts")
      expect_false(is.null(dimnames(submodel_data[[i]][["endogen"]])[[2]]))
      # Weight series must cover every sub-model of the data set.
      expect_setequal(dimnames(submodel_data[[i]][["weights"]])[[2]],
                      names(submodel_data))
    }
  }
})

test_that("shipped data sets pass the internal consistency check", {
  for (name in c("gvar2019", "gvar2023", "dees2007")) {
    utils::data(list = name, package = "bgvars", envir = environment())
    submodel_data <- get(name, envir = environment())[["submodel_data"]]
    expect_silent(bgvars:::.check_submodeldata(submodel_data))
  }
})

test_that("a sub-model does not carry a weight on itself", {
  submodel_data <- gvar_data()
  for (i in names(submodel_data)) {
    expect_true(all(submodel_data[[i]][["weights"]][, i] == 0))
  }
})
