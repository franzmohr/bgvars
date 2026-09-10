# Posterior simulation of GVEC sub-models on a scaled error correction term.
#
# The series in the error correction term are on the scale of the data, and the
# prior on the cointegration space is isotropic, so it treats them as if they
# were comparably scaled. add_posterior_coefficients() therefore simulates on a
# scaled term and puts the draws back afterwards.

test_that("the scaling methods come from bvartools", {
  # bgvars used to carry its own copy of both, and the copy of the rescaling
  # method had drifted: it built the inverse with diag() of a matrix, which
  # extracts a diagonal rather than building one, and it transformed alpha,
  # which the transformation does not touch. A 'vecxsubmodel' inherits the
  # bvartools methods through its 'bvecmodel' class, so there is nothing to
  # keep in this package.
  expect_null(getS3method("scale_error_correction", "vecxsubmodel",
                          optional = TRUE))
  expect_null(getS3method("rescale_error_correction", "vecxsubmodel",
                          optional = TRUE))

  object <- gvec_estimated(iterations = 20, burnin = 10)
  expect_s3_class(object[["submodels"]][["US"]][[1]], "bvecmodel")
})

test_that("scaling puts the error correction term on a comparable footing", {
  object <- gvec_estimated(iterations = 20, burnin = 10)
  model <- object[["submodels"]][["US"]][[1]]

  scaled <- scale_error_correction(model)
  factors <- attr(scaled[["data"]][["train"]][["w"]], "scale")

  expect_false(is.null(factors))
  expect_equal(names(factors), dimnames(model[["data"]][["train"]][["w"]])[[2]])

  spread <- function(x) {
    sds <- apply(as.matrix(x), 2, stats::sd)
    max(sds) / min(sds)
  }
  expect_lt(spread(scaled[["data"]][["train"]][["w"]]),
            spread(model[["data"]][["train"]][["w"]]))

  # The trend is divided by its own standard deviation, so it ends up at one.
  expect_equal(stats::sd(scaled[["data"]][["train"]][["w"]][, "trend"]), 1)
})

test_that("the transformation leaves the error correction term unchanged", {
  # With D the diagonal matrix of scaling factors, alpha beta' w is the same as
  # alpha (D beta)' (D^-1 w). If that did not hold, a scaled simulation would
  # be of a different model.
  object <- gvec_estimated(iterations = 20, burnin = 10)
  model <- object[["submodels"]][["US"]][[1]]

  scaled <- scale_error_correction(model)
  factors <- diag(attr(scaled[["data"]][["train"]][["w"]], "scale"))

  k_ect <- ncol(model[["data"]][["train"]][["w"]])
  beta <- matrix(model[["posterior"]][["beta"]][["coeffs"]][1, ], k_ect)

  plain <- t(beta) %*% t(as.matrix(model[["data"]][["train"]][["w"]]))
  transformed <- t(factors %*% beta) %*%
    t(as.matrix(scaled[["data"]][["train"]][["w"]]))

  expect_equal(plain, transformed)
})

test_that("what comes back is on the scale of the data", {
  # Only the draws should differ between a scaled and an unscaled simulation.
  # The data of the model are the same either way.
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)

  scaled <- add_posterior_coefficients(object, scale_ect = TRUE)
  plain <- add_posterior_coefficients(object, scale_ect = FALSE)

  for (submodel in names(object[["submodels"]])) {

    w <- scaled[["submodels"]][[submodel]][[1]][["data"]][["train"]][["w"]]

    # The attribute is dropped by the rescaling, which is what stops a second
    # rescaling from transforming the model again.
    expect_null(attr(w, "scale"), info = submodel)

    expect_equal(w, plain[["submodels"]][[submodel]][[1]][["data"]][["train"]][["w"]],
                 info = submodel)
  }
})

test_that("the initial values are left on the scale of the data", {
  # add_initial_values() produces reduced rank maximum likelihood estimates
  # that are compared with published cointegrating vectors, so the simulation
  # must not leave a transformed copy behind.
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)

  before <- object[["submodels"]][["US"]][[1]][["initial"]][["beta"]]
  after <- add_posterior_coefficients(object)[["submodels"]][["US"]][[1]][["initial"]][["beta"]]

  expect_equal(after, before)
})

test_that("scaling can be turned off", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)

  result <- add_posterior_coefficients(object, scale_ect = FALSE)

  expect_s3_class(result[["submodels"]][["US"]], "modellist")
  expect_null(attr(result[["submodels"]][["US"]][[1]][["data"]][["train"]][["w"]],
                   "scale"))
  expect_equal(nrow(result[["submodels"]][["US"]][[1]][["posterior"]][["beta"]][["coeffs"]]),
               20L)

  expect_error(add_posterior_coefficients(object, scale_ect = "yes"),
               "must be a single logical value")
})

test_that("a scaled simulation still solves into a global model", {
  object <- gvec_to_gvar(gvec_estimated(iterations = 30, burnin = 10))

  gvec <- submodels_to_gvar(object)

  expect_s3_class(gvec, "bvarmodel")
  expect_equal(gvec[["model"]][["k"]], 6L)
  expect_false(any(is.na(gvec[["posterior"]][["a"]][["coeffs"]])))
})

test_that("a cointegration rank above one is transformed correctly", {
  # Regression test. The starting values of beta are stored stacked, as one
  # column of k_ect * r elements rather than as a k_ect by r matrix, so
  # applying the scaling factors to them without reshaping first happens to
  # work at rank one and is not conformable above it.
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 2,
                          iterations = 20, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)

  before <- object[["submodels"]][["US"]][[1]][["initial"]][["beta"]]
  result <- add_posterior_coefficients(object)
  model <- result[["submodels"]][["US"]][[1]]

  expect_equal(model[["model"]][["rank"]], 2L)
  expect_equal(model[["initial"]][["beta"]], before)

  k_ect <- ncol(model[["data"]][["train"]][["w"]])
  expect_equal(ncol(model[["posterior"]][["beta"]][["coeffs"]]), k_ect * 2L)
})
