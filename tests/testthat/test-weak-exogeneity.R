# Checks the weak exogeneity test of the foreign variables of a VEC sub-model:
# that its F statistic is the F statistic base R computes for the same
# restriction, that it detects a dependent variable built to respond to the
# error correction term and does not flag one built independently of it, that
# the sample and the degrees of freedom follow 'lags', and that it refuses the
# models for which the assumption has no content.
#
# The collection methods are checked for shape here and nothing else; what they
# map is the method for one sub-model, which is checked above them.

# The error correction terms of one draw, built from the sub-model's own data
# and not from the function under test.
ect_of <- function(model, draw = 1) {
  beta <- as.matrix(model[["posterior"]][["beta"]][["coeffs"]])
  unclass(model[["data"]][["train"]][["w"]]) %*%
    matrix(beta[draw, ], model[["model"]][["k_beta"]], model[["model"]][["rank"]])
}

# The marginal model's own regressors, built independently of the function:
# with p_endogen = 1 the conditional model carries no lagged domestic
# differences, so the controls are the lags of the foreign differences alone.
marginal_data <- function(model, lags = 1) {
  exogen <- model[["data"]][["original"]][["exogen"]]
  differences <- unclass(diff(exogen))
  rows <- seq.int(1 + lags, nrow(differences))
  controls <- do.call(cbind, lapply(seq_len(lags), function(l) {
    differences[rows - l, , drop = FALSE]
  }))
  list(dependent = differences[rows, , drop = FALSE], controls = controls,
       rows = rows)
}

test_that("the statistic is the F statistic of the same restriction", {

  set.seed(1)
  object <- gvec_estimated(iterations = 60, burnin = 20, r = 1)
  model <- object[["submodels"]][["US"]][[1]]

  # One draw, so that the median the function reports is that draw and can be
  # compared with a regression built by hand.
  beta <- as.matrix(model[["posterior"]][["beta"]][["coeffs"]])
  model[["posterior"]][["beta"]][["coeffs"]] <- beta[1, , drop = FALSE]

  result <- weak_exogeneity_test(model, lags = 1)

  pieces <- marginal_data(model, lags = 1)
  terms <- ect_of(model)[pieces[["rows"]], , drop = FALSE]

  for (i in seq_along(result[["variable"]])) {
    dependent <- pieces[["dependent"]][, i]
    restricted <- stats::lm(dependent ~ pieces[["controls"]])
    unrestricted <- stats::lm(dependent ~ terms + pieces[["controls"]])
    comparison <- stats::anova(restricted, unrestricted)
    expect_equal(result[["statistic"]][i], comparison[["F"]][2],
                 tolerance = 1e-8)
    expect_equal(result[["p_value"]][i], comparison[["Pr(>F)"]][2],
                 tolerance = 1e-8)
    expect_equal(result[["df2"]][i], unrestricted[["df.residual"]])
  }
})

test_that("a foreign variable that responds to the relation is detected", {

  set.seed(2)
  object <- gvec_estimated(iterations = 40, burnin = 10, r = 1)
  model <- object[["submodels"]][["US"]][[1]]

  # The error correction term of the first draw, and a sub-model whose foreign
  # variables are rebuilt so that the difference of the first follows that term
  # and the difference of the second does not. The levels inside 'w' are left
  # alone, so the object is no longer internally consistent; what is being
  # checked is the behaviour of the statistic, not of the model.
  terms <- as.vector(ect_of(model))
  exogen <- model[["data"]][["original"]][["exogen"]]
  periods <- nrow(exogen)
  responds <- c(0, 3 * terms + stats::rnorm(periods - 1, sd = 0.01))
  independent <- c(0, stats::rnorm(periods - 1, sd = 0.1))
  rebuilt <- exogen
  rebuilt[, 1] <- cumsum(responds)
  rebuilt[, 2] <- cumsum(independent)
  model[["data"]][["original"]][["exogen"]] <- rebuilt

  result <- weak_exogeneity_test(model, lags = 1)

  # The term was built from one draw of beta and the statistic is computed at
  # every draw, so the other draws see a slightly different relation and a few
  # of them need not reject. What is being checked is that the variable built to
  # respond is detected overwhelmingly and the one built independently is not.
  expect_gt(result[["statistic"]][1], 20)
  expect_lt(result[["p_value"]][1], 1e-4)
  expect_gt(result[["rejected"]][1], 0.9)

  expect_lt(result[["statistic"]][2], 10)
  expect_gt(result[["p_value"]][2], 0.01)
})

test_that("'lags' sets the sample and the degrees of freedom", {

  set.seed(3)
  object <- gvec_estimated(iterations = 30, burnin = 10, r = 1)
  model <- object[["submodels"]][["US"]][[1]]
  observations <- nrow(model[["data"]][["train"]][["y"]])

  none <- weak_exogeneity_test(model, lags = 0)
  one <- weak_exogeneity_test(model, lags = 1)
  two <- weak_exogeneity_test(model, lags = 2)

  expect_equal(unique(none[["observations"]]), observations)
  expect_equal(unique(one[["observations"]]), observations - 1)
  expect_equal(unique(two[["observations"]]), observations - 2)

  # Two more regressors per lag, since the sub-model has two foreign variables,
  # and one observation fewer, so the residual degrees of freedom fall by three.
  expect_equal(unique(two[["df2"]]), unique(one[["df2"]]) - 3)
  expect_equal(unique(one[["df1"]]), 1L)
  expect_equal(unique(one[["rank"]]), 1L)
})

test_that("the share of rejecting draws follows 'level'", {

  set.seed(4)
  object <- gvec_estimated(iterations = 100, burnin = 20, r = 1)
  model <- object[["submodels"]][["US"]][[1]]

  strict <- weak_exogeneity_test(model, level = 0.01)
  loose <- weak_exogeneity_test(model, level = 0.2)

  expect_true(all(loose[["rejected"]] >= strict[["rejected"]]))
  expect_equal(strict[["statistic"]], loose[["statistic"]])
})

test_that("the collection methods keep the shape of the models", {

  set.seed(5)
  object <- gvec_estimated(iterations = 30, burnin = 10, r = 1)

  results <- weak_exogeneity_test(object)
  expect_type(results, "list")
  expect_named(results, names(object[["submodels"]]))
  expect_length(results[["US"]], length(object[["submodels"]][["US"]]))
  expect_s3_class(results[["US"]][[1]], "data.frame")
  expect_equal(results[["US"]][[1]],
               weak_exogeneity_test(object[["submodels"]][["US"]][[1]]))

  from_list <- weak_exogeneity_test(object[["submodels"]][["US"]])
  expect_type(from_list, "list")
  expect_equal(from_list[[1]], results[["US"]][[1]])
})

test_that("the test refuses the models the assumption has no content for", {

  set.seed(6)

  rank_zero <- gvec_estimated(iterations = 30, burnin = 10, r = 0)
  expect_error(weak_exogeneity_test(rank_zero[["submodels"]][["US"]][[1]]),
               "rank zero")

  var_model <- gvar_estimated(iterations = 30, burnin = 10)
  expect_error(weak_exogeneity_test(var_model), "global VAR")
  expect_error(weak_exogeneity_test(var_model[["submodels"]][["US"]][[1]]),
               "no cointegrating relations")

  object <- gvec_estimated(iterations = 30, burnin = 10, r = 1)
  model <- object[["submodels"]][["US"]][[1]]

  unestimated <- model
  unestimated[["posterior"]] <- NULL
  expect_error(weak_exogeneity_test(unestimated),
               "add_posterior_coefficients")

  # One cointegration space per period, which a time varying specification
  # gives, is refused rather than reshaped.
  varying <- model
  varying[["posterior"]][["beta"]][["coeffs"]] <-
    cbind(model[["posterior"]][["beta"]][["coeffs"]],
          model[["posterior"]][["beta"]][["coeffs"]])
  expect_error(weak_exogeneity_test(varying), "one cointegration space per")

  expect_error(weak_exogeneity_test(model, lags = -1), "non-negative integer")
  expect_error(weak_exogeneity_test(model, lags = 1.5), "non-negative integer")
  expect_error(weak_exogeneity_test(model, level = 0), "between 0 and 1")
  expect_error(weak_exogeneity_test(model, level = 1), "between 0 and 1")
  expect_error(weak_exogeneity_test(model, lags = 1000),
               "more lags of the foreign differences")
})
