# Scaling the error correction terms of a GVEC model.
#
# The series in the error correction term are on the scale of the data, and the
# prior on the cointegration space is isotropic, so it treats them as if they
# were comparably scaled. Scaling is a deliberate step of the user's:
# scale_error_correction() before the simulation and rescale_error_correction()
# after it. add_posterior_coefficients() never scales on its own.

gvec_prepared <- function(r = 1, tvp = FALSE, initial = TRUE) {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted",
                          trend = if (tvp) NULL else "restricted",
                          r = r, tvp = tvp,
                          iterations = 20, burnin = 10)
  if (tvp) {
    object <- align_model_obs(object)
  }
  object <- add_priors(object,
                       coef = if (tvp) {
                         list(v_i = 1, v_i_det = 1 / 10, shape = 3,
                              rate = 1e-8, rate_det = 1e-8)
                       } else {
                         list(v_i = 0)
                       },
                       coint = if (tvp) list(rho = 0.999) else list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  if (initial) {
    object <- add_initial_values(object)
  }
  object
}

w_of <- function(object, submodel = "US") {
  object[["submodels"]][[submodel]][[1]][["data"]][["train"]][["w"]]
}

initial_of <- function(object, name = "beta", submodel = "US") {
  object[["submodels"]][[submodel]][[1]][["initial"]][[name]]
}

test_that("the scaling methods for sub-models come from bvartools", {
  # bgvars used to carry its own copy of both, and the copy of the rescaling
  # method had drifted: it built the inverse with diag() of a matrix, which
  # extracts a diagonal rather than building one, and it transformed alpha,
  # which the transformation does not touch. A 'vecxsubmodel' inherits the
  # bvartools methods through its 'bvecmodel' class, so what this package adds
  # is only the step over all sub-models of a 'gvecmodel'.
  expect_null(getS3method("scale_error_correction", "vecxsubmodel",
                          optional = TRUE))
  expect_null(getS3method("rescale_error_correction", "vecxsubmodel",
                          optional = TRUE))
  expect_false(is.null(getS3method("scale_error_correction", "gvecmodel",
                                   optional = TRUE)))
  expect_false(is.null(getS3method("rescale_error_correction", "gvecmodel",
                                   optional = TRUE)))

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

test_that("posterior simulation does not scale on its own", {
  object <- gvec_prepared()

  result <- add_posterior_coefficients(object)

  for (submodel in names(object[["submodels"]])) {
    expect_null(attr(w_of(result, submodel), "scale"), info = submodel)
    expect_equal(w_of(result, submodel), w_of(object, submodel), info = submodel)
  }
  expect_equal(nrow(result[["submodels"]][["US"]][[1]][["posterior"]][["beta"]][["coeffs"]]),
               20L)

  # A script written for the implicit scaling must not run unscaled unnoticed.
  expect_error(add_posterior_coefficients(object, scale_ect = TRUE),
               "no longer supported")
})

test_that("every sub-model is scaled, and rescaling undoes it", {
  object <- gvec_prepared()

  scaled <- scale_error_correction(object)
  expect_s3_class(scaled, "gvecmodel")

  for (submodel in names(object[["submodels"]])) {
    expect_s3_class(scaled[["submodels"]][[submodel]], "modellist")
    expect_false(is.null(attr(w_of(scaled, submodel), "scale")), info = submodel)
  }

  rescaled <- rescale_error_correction(scaled)
  for (submodel in names(object[["submodels"]])) {
    expect_null(attr(w_of(rescaled, submodel), "scale"), info = submodel)
    expect_equal(unclass(w_of(rescaled, submodel)), unclass(w_of(object, submodel)),
                 info = submodel)
  }

  # A model that was never scaled has nothing to rescale.
  expect_error(rescale_error_correction(object), "does not have an attribute 'scale'")
})

test_that("a scaled simulation comes back on the scale of the data and solves", {
  object <- gvec_prepared()

  result <- rescale_error_correction(add_posterior_coefficients(scale_error_correction(object)))

  for (submodel in names(object[["submodels"]])) {
    expect_null(attr(w_of(result, submodel), "scale"), info = submodel)
    expect_equal(unclass(w_of(result, submodel)), unclass(w_of(object, submodel)),
                 info = submodel)
  }

  gvar <- submodels_to_gvar(result)
  expect_s3_class(gvar, "bvarmodel")
  expect_false(any(is.na(gvar[["posterior"]][["a"]][["coeffs"]])))
})

test_that("starting values follow the scaling, whichever comes first", {
  # Initial values added before scaling are transformed with the data and put
  # back afterwards.
  object <- gvec_prepared()
  before <- initial_of(object)

  scaled <- scale_error_correction(object)
  factors <- attr(w_of(scaled), "scale")
  k_ect <- length(factors)
  expect_equal(as.numeric(initial_of(scaled)),
               as.numeric(diag(factors) %*% matrix(before, k_ect)))

  back <- rescale_error_correction(scaled)
  expect_equal(initial_of(back), before)
  expect_equal(dim(initial_of(back)), dim(before))

  # Initial values added after scaling are estimated on the scaled data, and
  # rescaling puts them on the scale of the data. The cointegration space does
  # not depend on the scale of the variables, so the vector is the same up to
  # its normalisation.
  late <- add_initial_values(scale_error_correction(gvec_prepared(initial = FALSE)))
  late <- rescale_error_correction(late)
  a <- as.numeric(initial_of(late))
  b <- as.numeric(before)
  expect_equal(abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2)), 1, tolerance = 1e-6)
})

test_that("the path and the state before the sample of time varying vectors follow as well", {
  object <- gvec_prepared(tvp = TRUE)
  path <- initial_of(object, "beta")
  state <- initial_of(object, "beta_init")
  expect_false(is.null(path))
  expect_false(is.null(state))

  scaled <- scale_error_correction(object)
  factors <- attr(w_of(scaled), "scale")
  k_ect <- length(factors)

  expect_equal(as.numeric(initial_of(scaled, "beta_init")),
               as.numeric(diag(factors) %*% matrix(state, k_ect)))
  expect_equal(as.numeric(initial_of(scaled, "beta")),
               as.numeric(diag(factors) %*% matrix(path, k_ect)))

  back <- rescale_error_correction(scaled)
  expect_equal(initial_of(back, "beta"), path)
  expect_equal(initial_of(back, "beta_init"), state)
})

test_that("a cointegration rank above one is transformed correctly", {
  # Regression test. The starting values of beta are stored stacked, as one
  # column of k_ect * r elements rather than as a k_ect by r matrix, so
  # applying the scaling factors to them without reshaping first happens to
  # work at rank one and is not conformable above it.
  object <- gvec_prepared(r = 2)
  before <- initial_of(object)

  result <- rescale_error_correction(add_posterior_coefficients(scale_error_correction(object)))
  model <- result[["submodels"]][["US"]][[1]]

  expect_equal(model[["model"]][["rank"]], 2L)
  expect_equal(model[["initial"]][["beta"]], before)

  k_ect <- ncol(model[["data"]][["train"]][["w"]])
  expect_equal(ncol(model[["posterior"]][["beta"]][["coeffs"]]), k_ect * 2L)
})

test_that("centring passes on to the sub-models and back", {
  # Centring alone leaves beta in the units of the series, so its starting
  # values stay as they are, while the constant takes up the means.
  for (tvp in c(FALSE, TRUE)) {
    object <- gvec_prepared(tvp = tvp)
    centred <- scale_error_correction(object, scale = FALSE, centre = TRUE)
    back <- rescale_error_correction(centred)

    for (submodel in names(object[["submodels"]])) {
      info <- paste0(submodel, ", tvp: ", tvp)
      expect_false(is.null(attr(w_of(centred, submodel), "centre")), info = info)
      expect_null(attr(w_of(centred, submodel), "scale"), info = info)
      expect_equal(initial_of(centred, "beta", submodel), initial_of(object, "beta", submodel),
                   info = info)

      expect_null(attr(w_of(back, submodel), "centre"), info = info)
      expect_equal(unclass(w_of(back, submodel)), unclass(w_of(object, submodel)),
                   ignore_attr = TRUE, info = info)
      expect_equal(initial_of(back, "a", submodel), initial_of(object, "a", submodel), info = info)
      expect_equal(initial_of(back, "beta", submodel), initial_of(object, "beta", submodel),
                   info = info)
    }
  }
})

test_that("a centred and scaled GVEC model can be simulated and solved", {
  object <- gvec_prepared()
  result <- rescale_error_correction(add_posterior_coefficients(
    scale_error_correction(object, centre = TRUE)))

  for (submodel in names(object[["submodels"]])) {
    expect_null(attr(w_of(result, submodel), "centre"), info = submodel)
    expect_null(attr(w_of(result, submodel), "scale"), info = submodel)
    expect_equal(initial_of(result, "beta", submodel), initial_of(object, "beta", submodel),
                 info = submodel)
  }

  gvar <- submodels_to_gvar(result)
  expect_s3_class(gvar, "bvarmodel")
  expect_false(any(is.na(gvar[["posterior"]][["a"]][["coeffs"]])))
})
