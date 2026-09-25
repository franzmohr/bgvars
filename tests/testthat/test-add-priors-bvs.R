# BVS priors when the coefficient prior leaves out v_i or v_i_det.
#
# The warning about a flat prior compared coef$v_i and coef$v_i_det with zero
# directly. A Minnesota prior needs neither, and a VECX sub-model never filled
# in v_i_det, so both stopped with "argument is of length zero" before any
# prior was built.

bvs_kappa <- list(kappa1 = 2, kappa2 = 0.5, kappa3 = 0.5, kappa4 = 5)

bvs_varx <- function() {
  create_varxsubmodel(gvar_object(), "US",
                      endogen = c("y", "Dp"), p_endogen = 1,
                      exogen = c("y", "Dp"), p_exogen = 1,
                      varsel = "bvs", iterations = 10, burnin = 10)[[1]]
}

bvs_vecx <- function() {
  create_vecxsubmodel(gvec_object(), "US",
                      endogen = c("y", "Dp"), p_endogen = 2,
                      exogen = c("y", "Dp"), p_exogen = 2,
                      const = "unrestricted", r = 1,
                      varsel = "bvs", iterations = 10, burnin = 10)[[1]]
}

bvs_coint <- list(v_i = 0, p_tau_i = 1)
bvs_sigma <- list(df = 3, scale = 1)

test_that("BVS takes a Minnesota prior without v_i for a VARX sub-model", {
  model <- bvs_varx()
  expect_no_warning(
    prior <- add_priors(model, coef = list(minnesota = bvs_kappa),
                        sigma = bvs_sigma, varsel = list(inprior = 0.5))
  )
  minn <- minnesota_prior(model, kappa1 = 2, kappa2 = 0.5, kappa3 = 0.5,
                          kappa4 = 5, max_var = NULL, coint_var = FALSE,
                          sigma = "AR")
  expect_equal(prior[["priors"]][["a"]][["v_inv"]], minn[["v_inv"]])
  expect_false(is.null(prior[["priors"]][["a"]][["inprior"]]))
})

test_that("BVS takes a Minnesota prior without v_i for a VECX sub-model", {
  model <- bvs_vecx()
  expect_no_warning(
    prior <- add_priors(model, coef = list(minnesota = bvs_kappa),
                        coint = bvs_coint, sigma = bvs_sigma,
                        varsel = list(inprior = 0.5))
  )
  n_a <- ncol(model[["data"]][["train"]][["z"]])
  expect_equal(dim(prior[["priors"]][["a"]][["v_inv"]]), c(n_a, n_a))
  expect_true(all(diag(prior[["priors"]][["a"]][["v_inv"]]) > 0))
  expect_false(is.null(prior[["priors"]][["a"]][["inprior"]]))
})

test_that("BVS on a VECX sub-model uses v_i for the deterministic terms when v_i_det is not given", {
  model <- bvs_vecx()
  expect_no_warning(
    prior <- add_priors(model, coef = list(v_i = 2), coint = bvs_coint,
                        sigma = bvs_sigma, varsel = list(inprior = 0.5))
  )
  n_a <- ncol(model[["data"]][["train"]][["z"]])
  expect_equal(prior[["priors"]][["a"]][["v_inv"]], diag(2, n_a))
})

test_that("BVS with a flat prior warns, unless only flat deterministic terms are excluded", {
  expect_warning(
    add_priors(bvs_vecx(), coef = list(v_i = 0), coint = bvs_coint,
               sigma = bvs_sigma, varsel = list(inprior = 0.5)),
    "uninformative prior"
  )
  expect_warning(
    add_priors(bvs_varx(), coef = list(v_i = 1, v_i_det = 0),
               sigma = bvs_sigma, varsel = list(inprior = 0.5)),
    "uninformative prior"
  )
  expect_no_warning(
    add_priors(bvs_varx(), coef = list(v_i = 1, v_i_det = 0),
               sigma = bvs_sigma,
               varsel = list(inprior = 0.5, exclude_det = TRUE))
  )
})
