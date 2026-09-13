# What argument 'period' of add_weight_matrices() means.

test_that("a single period of the weight series gives constant weights of that period", {
  # A single number used to be read as the length of a rolling window, so
  # period = 2001 asked for a window of 2001 years and failed with a missing
  # value in a condition.
  submodel_data <- gvec_data()
  weights <- bgvars:::.create_weights_submodel(submodel_data[["US"]], period = 2001)

  trade <- submodel_data[["US"]][["weights"]]
  year <- as.numeric(trade[stats::time(trade) == 2001, ])
  names(year) <- colnames(trade)
  expect_equal(as.numeric(weights[1, ]), unname(year / sum(year)))
  expect_equal(nrow(unique(as.matrix(weights))), 1L)

  object <- create_gvecmodel(submodel_data = submodel_data)
  object <- add_weight_matrices(object, submodel_data, period = 2001)
  w <- get_weight_matrix(object, "US")
  n_endogen <- sum(object[["global"]][["index"]][, "submodel"] == "US")
  # Weakly exogenous output of the US: Japan and Canada in proportion to their
  # trade with the US in 2001.
  expect_equal(unname(w[n_endogen + 1, c("JP_y", "CA_y")]),
               unname(year[c("JP", "CA")] / sum(year[c("JP", "CA")])))

  # Short rolling windows keep their meaning.
  expect_gt(nrow(unique(as.matrix(
    bgvars:::.create_weights_submodel(submodel_data[["US"]], period = 3)))), 1L)
})

test_that("periods missing from the weight series and overlong windows are refused", {
  # Constant weights over periods the series does not have used to be summed
  # over the periods it does have, without a word, or became NaN.
  submodel_data <- gvec_data()

  expect_error(bgvars:::.create_weights_submodel(submodel_data[["US"]],
                                                 period = 2002:2004, name = "US"),
               "sub-model US runs from 1980 to 2003 and does not contain the period\\(s\\) 2004")
  expect_error(bgvars:::.create_weights_submodel(submodel_data[["US"]],
                                                 period = 30, name = "US"),
               "rolling window of 30 periods")

  object <- create_gvecmodel(submodel_data = submodel_data)
  expect_error(add_weight_matrices(object, submodel_data, period = 1978:1980),
               "sub-model US")
})
