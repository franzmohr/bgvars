# The prior on the cointegration space of a sub-model.
#
# add_priors() for a sub-model hands 'coint' to bvartools::cointspace_prior(),
# which builds the prior for the VEC models of bvartools as well: numeric,
# full-matrix and ML-centred locations, the loading shrinkage from the ML
# loadings, and the ML-centred transition of time varying cointegration. What
# has to hold for a sub-model is that its error correction term, which also
# holds the weakly exogenous and global variables, sets the dimensions. The
# sub-models below have lags of the weakly exogenous variables, a global
# variable and a restricted trend on purpose: that is where model$m counts
# differently from the error correction term.

vecx_ml_object <- function(tvp = FALSE, scale = FALSE) {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 3,
                          global = if (tvp) NULL else "poil",
                          s = if (tvp) NULL else 2,
                          const = "unrestricted", trend = "restricted", r = 1,
                          tvp = tvp, iterations = 20, burnin = 10)
  object <- align_model_obs(object)
  if (scale) {
    object <- scale_error_correction(object)
  }
  object
}

vecx_coef <- function(model) {
  if (model[["model"]][["tvp"]]) {
    list(v_i = 1, v_i_det = 0.1, shape = 3, rate = 1e-4)
  } else {
    list(v_i = 1, v_i_det = 0.1)
  }
}

# Priors for a 'gvecmodel' or a single 'vecxsubmodel'. Checks that warn are run
# on a single sub-model, so that the warning is raised once.
vecx_ml_priors <- function(object, coint) {
  model <- if (inherits(object, "gvecmodel")) object[["submodels"]][[1]][[1]] else object
  add_priors(object, coef = vecx_coef(model), coint = coint, sigma = list(df = 3, scale = 0.0001))
}

first_submodel <- function(object) object[["submodels"]][[1]][[1]]

submodel_list <- function(object) lapply(object[["submodels"]], function(x) x[[1]])

# The orthonormal basis of the space spanned by Johansen's estimate and of its
# orthogonal complement.
ml_space <- function(model) {
  beta <- bvartools:::.coint_ml(model)[["beta"]]
  basis <- qr.Q(qr(beta), complete = TRUE)
  r <- ncol(beta)
  list(h = basis[, seq_len(r), drop = FALSE],
       h_perp = basis[, -seq_len(r), drop = FALSE])
}

test_that("the prior is the one bvartools builds for the sub-model", {
  object <- vecx_ml_object(scale = TRUE)
  coints <- list(list(v_i = 0, p_tau_i = 1),
                 list(v_i = "ml", p_tau_i = "ml"),
                 list(v_i = 0.01, p_tau_i = "ml", weight = 10))
  for (coint in coints) {
    priors <- vecx_ml_priors(object, coint)
    for (i in seq_along(object[["submodels"]])) {
      model <- object[["submodels"]][[i]][[1]]
      expect_identical(priors[["submodels"]][[i]][[1]][["priors"]][["beta"]],
                       bvartools::cointspace_prior(model, coint))
    }
  }
})

test_that("a sub-model's prior is sized from its error correction term", {
  object <- vecx_ml_object()
  for (model in submodel_list(object)) {
    k_w <- ncol(model[["data"]][["train"]][["w"]])
    # The regressors model$m counts are not the series of the error correction
    # term, which is what used to give a numeric p_tau_i the wrong size.
    expect_false(model[["model"]][["m"]] + model[["model"]][["k"]] +
                   model[["model"]][["n_restricted"]] == k_w)
    prior <- vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = 1))[["priors"]][["beta"]]
    expect_identical(dim(prior[["p_tau_inv"]]), c(k_w, k_w))
  }
})

test_that("a cointegration space prior can be centred on the ML estimate", {
  object <- vecx_ml_priors(vecx_ml_object(scale = TRUE), list(v_i = 0.01, p_tau_i = "ml"))

  for (model in submodel_list(object)) {
    p_tau_inv <- model[["priors"]][["beta"]][["p_tau_inv"]]
    space <- ml_space(model)
    k_w <- ncol(model[["data"]][["train"]][["w"]])

    expect_equal(model[["priors"]][["beta"]][["v_inv"]], 0.01)
    expect_identical(dim(p_tau_inv), c(k_w, k_w))
    expect_equal(p_tau_inv, t(p_tau_inv))
    # Along the estimated space the prior is that of the uniform one ...
    expect_equal(p_tau_inv %*% space[["h"]], space[["h"]])
    # ... and away from it at least as tight.
    expect_true(all(eigen(p_tau_inv, symmetric = TRUE)$values >= 1 - 1e-8))
  }
})

test_that("the weight of the ML prior scales its precision off the space", {
  object <- vecx_ml_object(scale = TRUE)
  off_space <- function(weight, i) {
    model <- vecx_ml_priors(object, list(v_i = 0.01, p_tau_i = "ml", weight = weight))[["submodels"]][[i]][[1]]
    space <- ml_space(model)
    crossprod(space[["h_perp"]], model[["priors"]][["beta"]][["p_tau_inv"]] %*% space[["h_perp"]])
  }
  for (i in seq_along(object[["submodels"]])) {
    # Directions whose precision is capped at the uniform prior do not scale,
    # so the comparison is made where the prior is informative even at weight 1.
    strong <- eigen(off_space(10, i), symmetric = TRUE)$values
    weak <- eigen(off_space(1, i), symmetric = TRUE)$values
    informative <- weak > 1 + 1e-8
    expect_true(any(informative))
    expect_equal(strong[informative], 10 * weak[informative])
  }
})

test_that("the loading shrinkage can be set from the ML loadings", {
  object <- vecx_ml_priors(vecx_ml_object(scale = TRUE), list(v_i = "ml", p_tau_i = 1))
  for (model in submodel_list(object)) {
    v_inv <- model[["priors"]][["beta"]][["v_inv"]]
    k_w <- ncol(model[["data"]][["train"]][["w"]])
    expect_length(v_inv, 1)
    expect_true(is.numeric(v_inv) && is.finite(v_inv) && v_inv > 0)
    expect_equal(model[["priors"]][["beta"]][["p_tau_inv"]], diag(1, k_w))
  }
})

test_that("a full p_tau_i matrix is taken as it is", {
  model <- first_submodel(vecx_ml_object())
  k_w <- ncol(model[["data"]][["train"]][["w"]])
  p_tau_i <- diag(seq_len(k_w))
  p_tau_i[1, 2] <- p_tau_i[2, 1] <- 0.5
  prior <- vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = p_tau_i))[["priors"]][["beta"]]
  expect_equal(prior[["p_tau_inv"]], p_tau_i)
})

test_that("a cointegration prior of a sub-model is checked", {
  model <- first_submodel(vecx_ml_object())
  k_w <- ncol(model[["data"]][["train"]][["w"]])
  # With zero shrinkage the sampler never sees p_tau_i.
  expect_error(vecx_ml_priors(model, list(v_i = 0, p_tau_i = "ml")), "coint[$]v_i")
  expect_error(vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = "mle")), "coint[$]p_tau_i")
  expect_error(vecx_ml_priors(model, list(v_i = "mle", p_tau_i = 1)), "coint[$]v_i")
  expect_error(vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = "ml", weight = 0)), "coint[$]weight")
  expect_warning(vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = 1, weight = 2)), "coint[$]weight")
  expect_error(vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = diag(k_w + 1))),
               paste0(k_w, " x ", k_w))
  asymmetric <- diag(k_w)
  asymmetric[1, 2] <- 1
  expect_error(vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = asymmetric)), "symmetric")
  # Elements that are not part of the prior are refused rather than ignored.
  expect_error(vecx_ml_priors(model, list(v_i = 0.1, p_tau_i = 1, tau = 1)), "not recognised")
})

test_that("scaling is refused once the space prior has a direction", {
  object <- vecx_ml_priors(vecx_ml_object(), list(v_i = 0.01, p_tau_i = "ml"))
  expect_error(scale_error_correction(object), "before 'add_priors'")
  # Scaling first and centring the prior afterwards is the intended order.
  expect_no_error(vecx_ml_priors(vecx_ml_object(scale = TRUE), list(v_i = 0.01, p_tau_i = "ml")))
})

test_that("a global model with ML-centred priors can be estimated", {
  object <- vecx_ml_priors(vecx_ml_object(scale = TRUE), list(v_i = "ml", p_tau_i = "ml"))
  object <- add_initial_values(object)
  set.seed(271828)
  object <- add_posterior_coefficients(object)
  for (model in submodel_list(object)) {
    beta <- as.matrix(model[["posterior"]][["beta"]][["coeffs"]])
    expect_equal(nrow(beta), 20)
    expect_true(all(is.finite(beta)))
  }
})

# --- time varying cointegration space centred on the ML estimate -------------
#
# On these short samples a weight of 0.1 already asks for more than rho = 0.999
# allows in some directions of some sub-models, so the floor warning is
# expected wherever it is not the subject of the test.

test_that("a time varying cointegration space of a sub-model can be centred on the ML estimate", {
  rho <- 0.999
  coint <- list(rho = rho, p_tau_i = "ml", weight = 0.1)
  object <- vecx_ml_object(tvp = TRUE, scale = TRUE)
  priors <- suppressWarnings(vecx_ml_priors(object, coint))

  for (i in seq_along(object[["submodels"]])) {
    model <- priors[["submodels"]][[i]][[1]]
    prior <- model[["priors"]][["beta"]]
    p_tau <- prior[["p_tau"]]
    space <- ml_space(model)
    h <- space[["h"]]
    k_w <- ncol(model[["data"]][["train"]][["w"]])

    expect_identical(prior, suppressWarnings(bvartools::cointspace_prior(object[["submodels"]][[i]][[1]], coint)))
    expect_identical(dim(p_tau), c(k_w, k_w))
    expect_equal(length(prior[["mu"]]), nrow(prior[["v_inv"]]))
    # Along the estimated space the transition is rho alone, and off it no more.
    expect_equal(p_tau %*% h, h)
    off <- eigen(crossprod(space[["h_perp"]], p_tau %*% space[["h_perp"]]), symmetric = TRUE)$values
    expect_true(all(off <= 1 + 1e-10))
    expect_true(any(off < 1 - 1e-8))
  }

  # The loadings keep the shrinkage bgvars gives them for time varying
  # cointegration, which is separate from the prior on the space.
  loadings <- diag(first_submodel(priors)[["priors"]][["a"]][["v_inv"]])[1:2]
  expect_equal(loadings, rep(1 / (1 - rho^2), 2))
})

test_that("a time varying ML prior of a sub-model is checked", {
  model <- first_submodel(vecx_ml_object(tvp = TRUE))
  expect_warning(vecx_ml_priors(model, list(rho = 0.999, p_tau_i = 1)), "coint[$]p_tau_i")
  expect_warning(vecx_ml_priors(model, list(rho = 0.999, weight = 2)), "coint[$]weight")
  expect_error(vecx_ml_priors(model, list(rho = 0.999, p_tau_i = "ml", weight = -1)), "coint[$]weight")

  # The transition has a direction, so scaling afterwards is refused too.
  centred <- suppressWarnings(vecx_ml_priors(vecx_ml_object(tvp = TRUE),
                                             list(rho = 0.999, p_tau_i = "ml", weight = 0.1)))
  expect_error(scale_error_correction(centred), "before 'add_priors'")
})

# G, the matrix the loadings' prior is scaled by, is the error covariance in
# every model but the constant VEC with stochastic volatility, where it is fixed
# for the run. 'coint$g_i' gives its inverse; bvartools builds it and BayesTS
# reads it from /priors/beta/g_inv, so what has to hold for a sub-model is that
# it arrives there, sized by the endogenous variables rather than by the error
# correction term.
vecx_sv_object <- function(tvp = FALSE) {
  object <- add_submodels(gvec_object(),
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", r = 1,
                          tvp = tvp, error = "sv", iterations = 20, burnin = 10)
  align_model_obs(object)
}

vecx_sv_priors <- function(object, coint) {
  coef <- list(v_i = 1, v_i_det = 0.1)
  if (isTRUE(first_submodel(object)[["model"]][["tvp"]])) {
    coef <- c(coef, list(shape = 3, rate = 1e-4))
  }
  add_priors(object, coef = coef, coint = coint,
             sigma = list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                          state_variance = 0.05, offset = 1e-8))
}

test_that("coint$g_i sets the G of a sub-model with stochastic volatility", {
  object <- vecx_sv_object()
  for (g_i in list(2, c(2, 3), "ml")) {
    coint <- list(v_i = 0, p_tau_i = 1, g_i = g_i)
    priors <- vecx_sv_priors(object, coint)
    for (i in seq_along(object[["submodels"]])) {
      model <- object[["submodels"]][[i]][[1]]
      g_inv <- priors[["submodels"]][[i]][[1]][["priors"]][["beta"]][["g_inv"]]
      expect_identical(priors[["submodels"]][[i]][[1]][["priors"]][["beta"]],
                       bvartools::cointspace_prior(model, coint))
      expect_identical(dim(g_inv), rep(model[["model"]][["k"]], 2))
    }
  }

  # Left out, G comes from the starting values and nothing is stored.
  priors <- vecx_sv_priors(object, list(v_i = 0, p_tau_i = 1))
  expect_null(first_submodel(priors)[["priors"]][["beta"]][["g_inv"]])
})

test_that("coint$g_i is refused where G is the error covariance", {
  expect_error(vecx_sv_priors(vecx_sv_object(tvp = TRUE), list(rho = 0.999, g_i = "ml")),
               "coint[$]g_i")
  expect_error(vecx_ml_priors(vecx_ml_object(), list(v_i = 0, p_tau_i = 1, g_i = 2)),
               "coint[$]g_i")
})

test_that("coint$g_i is written to the sub-model's file", {
  skip_if_not_installed("hdf5r")
  object <- vecx_sv_priors(vecx_sv_object(), list(v_i = 0, p_tau_i = 1, g_i = "ml"))
  object <- add_initial_values(object)
  folder <- tempfile("g-i-export")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder)

  file <- file.path(folder, "submodels", names(object[["submodels"]])[1], "001.h5")
  h5 <- hdf5r::H5File$new(file, mode = "r")
  on.exit(h5$close_all(), add = TRUE, after = FALSE)
  expect_equal(as.matrix(h5[["priors/beta/g_inv"]]$read()),
               first_submodel(object)[["priors"]][["beta"]][["g_inv"]],
               ignore_attr = TRUE)
})
