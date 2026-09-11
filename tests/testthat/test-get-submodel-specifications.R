test_that("get_submodel_specifications summarises the estimated sub-models", {
  object <- gvar_estimated(iterations = 20, burnin = 10, p_endogen = 1:2)

  specs <- get_submodel_specifications(object[["submodels"]][["US"]])

  expect_s3_class(specs, "data.frame")
  expect_equal(nrow(specs), 2L)
  expect_true(all(specs[, "type"] == "VARX"))
  expect_equal(specs[, "lag_domestic"], 1:2)

  expect_true(all(specs[, "var_domestic"] == "y, Dp"))

  # The models of one sub-model carry no names, so there is no group to report
  # and the column is dropped.
  expect_false("group" %in% names(specs))
})

test_that("the type reports what makes a sub-model non-standard", {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 0,
                          global = "poil", s = 0,
                          tvp = TRUE, error = "sv+covar",
                          iterations = 10, burnin = 5)
  object <- align_model_obs(object)

  specs <- get_submodel_specifications(object[["submodels"]][["US"]])

  # A time varying model with stochastic volatility says so in one string.
  expect_equal(specs[, "type"], "TVP-SV-VARX")
  expect_equal(specs[, "var_global"], "poil")
})

test_that("columns that apply to no model are dropped", {
  object <- gvar_estimated(iterations = 20, burnin = 10)

  specs <- get_submodel_specifications(object[["submodels"]][["US"]])

  # A VARX sub-model has no cointegration rank, so the column goes.
  expect_false("r" %in% names(specs))
})
