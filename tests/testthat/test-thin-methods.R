# Thinning GVAR and GVEC models, and what the global model keeps of it.
#
# bvartools labels the draws it keeps by where they sit in the chain, and every
# object it derives from them keeps those labels. thin() for the objects of this
# package forwards to the sub-models, and submodels_to_gvar() stacks draw i of
# every sub-model into draw i of the global model, so the labels have to agree
# across sub-models and have to reach the global model.

test_that("thin thins every sub-model of a GVAR model alike", {
  object <- gvar_estimated(iterations = 12, burnin = 5)
  thinned <- thin(object, thin = 3)
  keep <- seq(3, 12, 3)

  expect_s3_class(thinned, "gvarmodel")
  for (submodel in names(object[["submodels"]])) {
    expect_s3_class(thinned[["submodels"]][[submodel]], "modellist")
    before <- object[["submodels"]][[submodel]][[1]][["posterior"]]
    after <- thinned[["submodels"]][[submodel]][[1]][["posterior"]]
    for (block in c("a", "u_sigma_inv")) {
      expect_identical(c(as.matrix(after[[block]][["coeffs"]])),
                       c(as.matrix(before[[block]][["coeffs"]])[keep, , drop = FALSE]),
                       label = paste(submodel, block))
      expect_equal(attr(after[[block]][["coeffs"]], "mcpar"), c(3, 12, 3),
                   label = paste(submodel, block))
    }
  }
})

test_that("thin thins every sub-model of a GVEC model alike", {
  object <- gvec_estimated(iterations = 12, burnin = 5, r = 1)
  thinned <- thin(object, thin = 4)

  expect_s3_class(thinned, "gvecmodel")
  for (submodel in names(object[["submodels"]])) {
    for (block in c("a", "beta", "u_sigma_inv")) {
      draws <- thinned[["submodels"]][[submodel]][[1]][["posterior"]][[block]][["coeffs"]]
      expect_identical(nrow(draws), 3L, label = paste(submodel, block))
      expect_equal(attr(draws, "mcpar"), c(4, 12, 4), label = paste(submodel, block))
    }
  }
})

test_that("thin refuses sub-models that have not been estimated", {
  object <- add_submodels(gvar_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          iterations = 10, burnin = 5)
  expect_error(thin(object, thin = 2), "has not been estimated")
})

test_that("the global model keeps the labels of sub-models thinned by their sampler", {
  gvar <- submodels_to_gvar(gvar_estimated(iterations = 4, burnin = 5, thin = 3))

  for (draws in list(gvar[["posterior"]][["a"]][["coeffs"]],
                     gvar[["posterior"]][["u_sigma_inv"]][["coeffs"]],
                     gvar[["g"]])) {
    expect_equal(attr(draws, "mcpar"), c(3, 12, 3))
  }
  expect_identical(gvar[["model"]][["thin"]], 3L)
  expect_identical(gvar[["model"]][["iterations"]], 4L)
})

test_that("the global model keeps the labels of sub-models thinned afterwards", {
  gvar <- submodels_to_gvar(thin(gvar_estimated(iterations = 12, burnin = 5), thin = 3))

  expect_equal(attr(gvar[["posterior"]][["a"]][["coeffs"]], "mcpar"), c(3, 12, 3))
  expect_identical(gvar[["model"]][["thin"]], 3L)
})

test_that("the global model of unthinned sub-models carries no thinning interval", {
  gvar <- submodels_to_gvar(gvar_estimated(iterations = 12, burnin = 5))

  expect_equal(attr(gvar[["posterior"]][["a"]][["coeffs"]], "mcpar"), c(1, 12, 1))
  expect_null(gvar[["model"]][["thin"]])
})

test_that("sub-models whose draws are labelled differently are refused", {
  # As many draws each, but one sub-model kept every third draw of a longer chain.
  object <- gvar_estimated(iterations = 4, burnin = 5)
  thinned <- gvar_estimated(iterations = 4, burnin = 5, thin = 3)
  object[["submodels"]][["JP"]] <- thinned[["submodels"]][["JP"]]

  expect_error(submodels_to_gvar(object), "labelled differently")
})

test_that("a time varying sub-model keeps its labels at the period it is solved at", {
  object <- add_submodels(gvar_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 1,
                          deterministic = "const",
                          tvp = TRUE, error = "sv",
                          iterations = 4, burnin = 5, thin = 3)
  object <- align_model_obs(object)
  object <- add_priors(object,
                       coef = list(v_i = 1, v_i_det = 1 / 10, shape = 3,
                                   rate = 1e-8, rate_det = 1e-8),
                       sigma = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                                    state_variance = 0.05, offset = 1e-8))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  gvar <- submodels_to_gvar(object)
  expect_equal(attr(gvar[["posterior"]][["a"]][["coeffs"]], "mcpar"), c(3, 12, 3))
  expect_equal(attr(gvar[["posterior"]][["u_sigma_inv"]][["coeffs"]], "mcpar"), c(3, 12, 3))
  expect_identical(gvar[["model"]][["thin"]], 3L)
})
