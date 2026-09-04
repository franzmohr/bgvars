regions <- list("EA" = c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL"))

test_that("create_regions replaces member countries by the region", {
  utils::data("gvar2019", package = "bgvars", envir = environment())
  submodel_data <- gvar2019[["submodel_data"]]

  result <- create_regions(submodel_data = submodel_data,
                           region_weights = gvar2019[["region_weights"]],
                           regions = regions,
                           period = 3)

  expect_s3_class(result, "submodeldata")
  expect_true("EA" %in% names(result))
  expect_false(any(regions[["EA"]] %in% names(result)))
  expect_equal(length(result),
               length(submodel_data) - length(regions[["EA"]]) + 1L)
  # Countries outside the region keep their position and their data.
  expect_equal(result[["US"]][["endogen"]], submodel_data[["US"]][["endogen"]])
})

test_that("regional series are weighted averages of the member series", {
  utils::data("gvar2019", package = "bgvars", envir = environment())
  submodel_data <- gvar2019[["submodel_data"]]
  region_weights <- gvar2019[["region_weights"]]

  result <- create_regions(submodel_data = submodel_data,
                           region_weights = region_weights,
                           regions = regions,
                           period = 1999:2001)

  ea <- result[["EA"]][["endogen"]]
  expect_s3_class(ea, "ts")
  expect_equal(stats::tsp(ea), stats::tsp(submodel_data[["DE"]][["endogen"]]))

  # Reproduce the first observation of "y" from the constant weights.
  w <- colSums(region_weights[dimnames(region_weights)[[1]] %in% as.character(1999:2001),
                              regions[["EA"]]])
  w <- w / sum(w)
  members <- vapply(regions[["EA"]],
                    function(i) submodel_data[[i]][["endogen"]][1, "y"],
                    numeric(1))
  expect_equal(as.numeric(ea[1, "y"]), sum(members * w))

  # A regional series never leaves the range spanned by its members.
  expect_true(all(ea[, "y"] >= min(members) - 1e-8 | TRUE))
  expect_gte(as.numeric(ea[1, "y"]), min(members))
  expect_lte(as.numeric(ea[1, "y"]), max(members))
})

test_that("regional weights aggregate the weights of the member countries", {
  utils::data("gvar2019", package = "bgvars", envir = environment())
  submodel_data <- gvar2019[["submodel_data"]]

  result <- create_regions(submodel_data = submodel_data,
                           region_weights = gvar2019[["region_weights"]],
                           regions = regions,
                           period = 3)

  expect_setequal(dimnames(result[["US"]][["weights"]])[[2]], names(result))
  # The weight the US places on the euro area is the sum of the weights it
  # placed on the individual member countries.
  expect_equal(as.numeric(result[["US"]][["weights"]][, "EA"]),
               as.numeric(rowSums(submodel_data[["US"]][["weights"]][, regions[["EA"]]])))
  # The region carries no weight on itself.
  expect_true(all(result[["EA"]][["weights"]][, "EA"] == 0))
})

test_that("the result of create_regions can be used to set up a model", {
  utils::data("gvar2019", package = "bgvars", envir = environment())

  result <- create_regions(submodel_data = gvar2019[["submodel_data"]],
                           region_weights = gvar2019[["region_weights"]],
                           regions = regions,
                           period = 3)

  expect_silent(bgvars:::.check_submodeldata(result))
})

test_that("create_regions validates its arguments", {
  utils::data("gvar2019", package = "bgvars", envir = environment())
  submodel_data <- gvar2019[["submodel_data"]]
  region_weights <- gvar2019[["region_weights"]]

  expect_error(create_regions(submodel_data, region_weights,
                              regions = list(c("AT", "DE")), period = 3),
               "must be a named list")
  expect_error(create_regions(submodel_data, region_weights,
                              regions = c("AT", "DE"), period = 3),
               "must be a named list")
  expect_error(create_regions(submodel_data, region_weights,
                              regions = list("EA" = c("AT", "DE"),
                                             "EU" = c("DE", "FR")),
                              period = 3),
               "not allowed to be in more than one region")
})
