# Time varying parameters and Bayesian variable selection for VECX sub-models.
#
# The samplers of bvartools have supported both for error correction models all
# along; what was missing here were the priors that feed them.

vecx_tvp <- function(varsel = "bvs", error = "sv+covar", tvp = TRUE, r = 1) {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", r = r,
                          tvp = tvp, error = error, varsel = varsel,
                          iterations = 20, burnin = 10)
  align_model_obs(object)
}

vecx_priors <- function(object, ...) {
  coint <- if (object[["submodels"]][[1]][[1]][["model"]][["tvp"]]) {
    list(rho = 0.999)
  } else {
    list(v_i = 0, p_tau_i = 1)
  }
  coef <- list(v_i = 1, v_i_det = 1 / 10, shape = 3, rate = 1e-8, rate_det = 1e-6)
  sigma <- list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
                state_variance = 0.05, offset = 1e-8)
  add_priors(object, coef = coef, coint = coint, sigma = sigma, ...)
}

test_that("the loadings are excluded from the inclusion prior", {
  model <- vecx_tvp()[["submodels"]][["US"]][[1]]
  prior <- inclusion_prior(model, prob = 0.5, exclude_deterministics = TRUE)

  k <- model[["model"]][["k_endogen"]]
  r <- model[["model"]][["rank"]]
  n_alpha <- k * r
  n_det <- k * model[["model"]][["n"]]
  n_z <- ncol(model[["data"]][["train"]][["z"]])

  # One prior per coefficient, and the loadings carry no inclusion probability
  # because the cointegration term is drawn in a step of its own.
  expect_length(prior[["prior"]], n_z)
  expect_true(all(is.na(prior[["prior"]][1:n_alpha])))
  expect_false(any(is.na(prior[["prior"]][-(1:n_alpha)])))

  # Neither the loadings nor the deterministic terms are selected over.
  expect_length(prior[["include"]], n_z - n_alpha - n_det)
  expect_equal(min(prior[["include"]]), n_alpha + 1)
})

test_that("deterministic terms stay in the selection when asked for", {
  model <- vecx_tvp()[["submodels"]][["US"]][[1]]
  prior <- inclusion_prior(model, prob = 0.5, exclude_deterministics = FALSE)

  k <- model[["model"]][["k_endogen"]]
  n_alpha <- k * model[["model"]][["rank"]]
  n_z <- ncol(model[["data"]][["train"]][["z"]])

  expect_length(prior[["include"]], n_z - n_alpha)
})

test_that("a rank of zero leaves every coefficient selectable", {
  model <- vecx_tvp(r = 0)[["submodels"]][["US"]][[1]]
  prior <- inclusion_prior(model, prob = 0.5, exclude_deterministics = FALSE)

  expect_false(any(is.na(prior[["prior"]])))
  expect_length(prior[["include"]], ncol(model[["data"]][["train"]][["z"]]))
})

test_that("Minnesota-like inclusion priors weight own lags differently", {
  model <- vecx_tvp()[["submodels"]][["US"]][[1]]
  prior <- inclusion_prior(model, exclude_deterministics = FALSE,
                           minnesota_like = TRUE,
                           kappa1 = 0.8, kappa2 = 0.5, kappa3 = 0.4, kappa4 = 0.9)

  k <- model[["model"]][["k_endogen"]]
  n_alpha <- k * model[["model"]][["rank"]]
  n_det <- k * model[["model"]][["n"]]
  values <- as.numeric(prior[["prior"]])

  # One lag of the differenced endogenous variables: own lags get kappa1 and
  # the rest kappa2, laid out as vec of a k x k block.
  endogen_block <- matrix(values[n_alpha + 1:(k * k)], k)
  expect_equal(diag(endogen_block), rep(0.8, k))
  expect_equal(endogen_block[upper.tri(endogen_block)], 0.5)

  # Deterministic terms get kappa4.
  expect_equal(tail(values, n_det), rep(0.9, n_det))
})

test_that("a time varying VECX sub-model gets a state equation for beta", {
  object <- vecx_priors(vecx_tvp(varsel = "none"))
  priors <- object[["submodels"]][["US"]][[1]][["priors"]]
  model <- object[["submodels"]][["US"]][[1]][["model"]]

  rho <- 0.999
  n_beta <- model[["rank"]] * model[["k_beta"]]

  expect_equal(priors[["beta"]][["type"]], "cointspace")
  expect_equal(priors[["beta"]][["rho"]], rho)
  expect_equal(dim(priors[["beta"]][["v_inv"]]), c(n_beta, n_beta))

  # The prior on the state before the sample is the stationary distribution of
  # beta_t = rho beta_{t-1} + eta_t, whose precision is (1 - rho^2) I.
  expect_equal(diag(priors[["beta"]][["v_inv"]]), rep(1 - rho^2, n_beta))

  # A time varying model has no cointegration space prior.
  expect_null(priors[["beta"]][["p_tau_inv"]])
})

test_that("the loadings compensate the scale of the cointegration space", {
  object <- vecx_priors(vecx_tvp(varsel = "none"))
  submodel <- object[["submodels"]][["US"]][[1]]
  priors <- submodel[["priors"]]
  model <- submodel[["model"]]

  rho <- 0.999
  n_alpha <- model[["k_endogen"]] * model[["rank"]]
  n_z <- ncol(submodel[["data"]][["train"]][["z"]])
  precision <- diag(priors[["a"]][["v_inv"]])

  # Only alpha beta' is identified, so the prior variance of the loadings is
  # shrunk by the factor that the stationary variance of beta inflates it.
  expect_equal(precision[1:n_alpha], rep(1 / (1 - rho^2), n_alpha))
  # The remaining non-deterministic coefficients keep coef$v_i.
  expect_equal(precision[n_alpha + 1], 1)

  # The state equation of the coefficients is given a prior of its own.
  expect_length(priors[["a"]][["shape"]], n_z)
  expect_length(priors[["a"]][["rate"]], n_z)
})

test_that("a time varying VECX sub-model can be estimated with BVS", {
  object <- vecx_priors(vecx_tvp(), varsel = list(inprior = 0.5, exclude_det = TRUE))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  model <- object[["submodels"]][["US"]][[1]]
  tt <- nrow(model[["data"]][["train"]][["y"]])
  n_z <- ncol(model[["data"]][["train"]][["z"]])
  n_alpha <- model[["model"]][["k_endogen"]] * model[["model"]][["rank"]]

  # The coefficients are a path, the inclusion parameters are not: one is drawn
  # per coefficient, so a regressor is in the sub-model for the whole sample.
  expect_equal(ncol(model[["posterior"]][["a"]][["coeffs"]]), n_z * tt)
  expect_equal(ncol(model[["posterior"]][["a"]][["lambda"]]), n_z)
  expect_equal(ncol(model[["posterior"]][["beta"]][["coeffs"]]),
               model[["model"]][["rank"]] * model[["model"]][["k_beta"]] * tt)

  # The loadings are never switched off.
  lambda <- colMeans(model[["posterior"]][["a"]][["lambda"]])
  expect_equal(lambda[1:n_alpha], rep(1, n_alpha))

  expect_no_error(summary(model))
})

test_that("variable selection over the coefficients alone needs time varying parameters", {
  object <- vecx_tvp(tvp = FALSE, error = "sv+covar")

  # A constant coefficient sampler reads one selection scheme for coefficients
  # and covariances together.
  expect_error(vecx_priors(object, varsel = list(inprior = 0.5, exclude_det = TRUE)),
               "applies one selection scheme to both")

  expect_no_error(vecx_priors(object, varsel = list(inprior = 0.5, exclude_det = TRUE,
                                                    covar = TRUE)))
})

test_that("SSVS is refused for error correction sub-models", {
  object <- vecx_tvp(varsel = "ssvs")

  expect_error(vecx_priors(object, varsel = list(inprior = 0.5, tau = c(0.1, 10))),
               "not implemented for error correction sub-models")
})
