# The two parameterisations of the state variances of a sub-model with time
# varying parameters: the gamma prior on their precisions that 'shape' and 'rate'
# ask for, and the normal prior on their standard deviations that 'omega_v' asks
# for. What is checked here is that each is stored where the samplers look for it,
# that the two are never stored together, and that a 'coef' naming both is
# refused. Whether BayesTS then draws them is not this file's business;
# test-gvec-pipeline.R and test-gvar-pipeline.R estimate.

tvp_vecx <- function(r = 1, error = "sv+covar") {
  suppressMessages(
    create_vecxsubmodel(gvec_object(), submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 2,
                        exogen = c("y", "Dp"), p_exogen = 2,
                        const = "unrestricted", r = r, tvp = TRUE, error = error,
                        iterations = 10, burnin = 10))[[1]]
}

tvp_varx <- function(error = "sv+covar") {
  suppressMessages(
    create_varxsubmodel(gvar_object(), submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 2,
                        exogen = c("y", "Dp"), p_exogen = 2,
                        tvp = TRUE, error = error,
                        iterations = 10, burnin = 10))[[1]]
}

sv_prior <- function() {
  list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
       state_variance = 0.05, offset = 1e-8)
}

tvp_coint <- function() {
  list(rho = 0.999, rho_min = 0.999, rho_max = 1)
}

test_that("shape and rate put a gamma prior on the state precisions", {
  result <- add_priors(tvp_vecx(),
                       coef = list(v_i = 1, shape = 3, rate = 1e-4),
                       coint = tvp_coint(), sigma = sv_prior())

  expect_false(is.null(result[["priors"]][["a"]][["shape"]]))
  expect_false(is.null(result[["priors"]][["a"]][["rate"]]))
  expect_null(result[["priors"]][["a"]][["omega_v"]])
  expect_null(result[["priors"]][["psi"]][["omega_v"]])
})

test_that("omega_v puts a normal prior on the state standard deviations", {
  model <- tvp_vecx()
  result <- add_priors(model,
                       coef = list(v_i = 1, omega_v = 0.01),
                       coint = tvp_coint(), sigma = sv_prior())

  a <- result[["priors"]][["a"]]
  expect_false(is.null(a[["omega_v"]]))
  # The two parameterisations are mutually exclusive, and BayesTS rejects a block
  # that carries both, so the one not asked for must be absent rather than left
  # at some default.
  expect_null(a[["shape"]])
  expect_null(a[["rate"]])
  expect_equal(as.numeric(a[["omega_v"]]), rep(0.01, nrow(a[["omega_v"]])))

  # The covariance block follows the coefficients, as it does for shape and rate.
  k <- model[["model"]][["k"]]
  psi <- result[["priors"]][["psi"]]
  expect_equal(nrow(psi[["omega_v"]]), k * (k - 1) / 2)
  expect_null(psi[["shape"]])
  expect_null(psi[["rate"]])
})

test_that("omega_v_alpha is the loadings' own prior", {
  r <- 1
  model <- tvp_vecx(r = r)
  result <- add_priors(model,
                       coef = list(v_i = 1, omega_v = 0.01, omega_v_alpha = 1e-8),
                       coint = tvp_coint(), sigma = sv_prior())

  omega_v <- as.numeric(result[["priors"]][["a"]][["omega_v"]])
  n_alpha <- model[["model"]][["k"]] * r
  expect_equal(omega_v[seq_len(n_alpha)], rep(1e-8, n_alpha))
  expect_equal(omega_v[-seq_len(n_alpha)], rep(0.01, length(omega_v) - n_alpha))
})

test_that("the two parameterisations of the state variances cannot be combined", {
  expect_error(
    add_priors(tvp_vecx(), coef = list(v_i = 1, shape = 3, rate = 1e-4, omega_v = 0.01),
               coint = tvp_coint(), sigma = sv_prior()),
    "gives both")
  # The message names the elements that were given beside omega_v, so that a
  # 'rate_det' left over from an earlier specification is findable.
  expect_error(
    add_priors(tvp_vecx(), coef = list(v_i = 1, omega_v = 0.01, rate_det = 1e-10),
               coint = tvp_coint(), sigma = sv_prior()),
    "rate_det")
})

test_that("omega_v must be positive and omega_v_alpha needs omega_v", {
  expect_error(
    add_priors(tvp_vecx(), coef = list(v_i = 1, omega_v = 0),
               coint = tvp_coint(), sigma = sv_prior()),
    "must be positive")
  expect_error(
    add_priors(tvp_vecx(), coef = list(v_i = 1, omega_v = 0.01, omega_v_alpha = -1),
               coint = tvp_coint(), sigma = sv_prior()),
    "omega_v_alpha' must be positive")
  expect_error(
    add_priors(tvp_vecx(), coef = list(v_i = 1, shape = 3, rate = 1e-4,
                                       omega_v_alpha = 1e-8),
               coint = tvp_coint(), sigma = sv_prior()),
    "needs")
})

test_that("a TVP model still needs one of the two parameterisations", {
  expect_error(
    add_priors(tvp_vecx(), coef = list(v_i = 1),
               coint = tvp_coint(), sigma = sv_prior()),
    "omega_v")
})

test_that("a VARX sub-model takes omega_v in both blocks", {
  model <- tvp_varx()
  result <- add_priors(model, coef = list(v_i = 1, omega_v = 0.01),
                       sigma = sv_prior())

  expect_false(is.null(result[["priors"]][["a"]][["omega_v"]]))
  expect_null(result[["priors"]][["a"]][["shape"]])
  expect_false(is.null(result[["priors"]][["psi"]][["omega_v"]]))
  expect_null(result[["priors"]][["psi"]][["shape"]])

  # Only an error correction sub-model has loadings, so the element that scales
  # them is not recognised here.
  expect_error(
    add_priors(tvp_varx(), coef = list(v_i = 1, omega_v = 0.01, omega_v_alpha = 1e-8),
               sigma = sv_prior()),
    "not recognised")
})
