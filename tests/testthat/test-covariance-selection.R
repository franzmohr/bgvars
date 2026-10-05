# Whether variable selection reaches the error covariance block is a choice of
# its own, made with `varsel$covar`, and it is available for every sub-model
# that has such a block rather than only for the time varying ones. Selecting
# over the coefficients and not over the covariances used to be refused, because
# a constant coefficient sampler read one scheme for the whole model.
#
# What the choice costs is checked in bvartools and in BayesTS, which draw the
# models; here it is the prior that is checked, and that both sub-model classes
# answer the argument the same way.

covar_submodels <- function(covar, tvp = FALSE, form = c("vec", "var")) {

  form <- match.arg(form)
  object <- if (form == "vec") gvec_object() else gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          tvp = tvp, error = "sv+covar", varsel = "bvs",
                          iterations = 10, burnin = 5,
                          r = if (form == "vec") 1 else NULL,
                          const = if (form == "vec") "unrestricted" else NULL,
                          trend = if (form == "vec") "restricted" else NULL,
                          deterministic = if (form == "var") "const" else NULL)

  coef <- list(v_i = 1, v_i_det = 1 / 10)
  if (tvp) {
    coef <- c(coef, list(shape = 3, rate = 0.0001))
  }
  coint <- if (form == "vec") {
    list(v_i = 1, p_tau_i = 1, rho = if (tvp) 0.999 else NULL)
  }
  args <- list(object, coef = coef,
               sigma = list(mu = 0, v_i = 1 / 100, shape = 3, rate = 0.0001,
                            state_variance = 0.05, offset = 1e-8),
               varsel = list(inprior = 0.5, exclude_det = TRUE, covar = covar))
  if (!is.null(coint)) {
    args <- c(args, list(coint = coint))
  }
  suppressWarnings(do.call(add_priors, args))
}

test_that("a sub-model may select over its coefficients and not its covariances", {

  for (form in c("vec", "var")) {

    selected <- covar_submodels(covar = TRUE, form = form)[["submodels"]][["US"]][[1]]
    expect_equal(selected[["priors"]][["psi"]][["varsel"]], "bvs")
    expect_false(is.null(selected[["priors"]][["psi"]][["inprior"]]))

    # This combination was refused until the samplers with constant
    # coefficients read the covariance block's own scheme.
    alone <- covar_submodels(covar = FALSE, form = form)[["submodels"]][["US"]][[1]]
    expect_equal(alone[["priors"]][["psi"]][["varsel"]], "none")
    expect_null(alone[["priors"]][["psi"]][["inprior"]])

    # The coefficients are selected over either way, which is the point of
    # being able to ask for one and not the other.
    expect_equal(alone[["model"]][["varsel"]], "bvs")
    expect_false(is.null(alone[["priors"]][["a"]][["inprior"]]))
  }
})

test_that("a time varying sub-model answers the argument the same way", {

  selected <- covar_submodels(covar = TRUE, tvp = TRUE)[["submodels"]][["US"]][[1]]
  alone <- covar_submodels(covar = FALSE, tvp = TRUE)[["submodels"]][["US"]][[1]]

  expect_equal(selected[["priors"]][["psi"]][["varsel"]], "bvs")
  expect_equal(alone[["priors"]][["psi"]][["varsel"]], "none")
})

test_that("one endogenous variable has no covariance block to select over", {

  object <- gvec_object()
  object <- add_submodels(object, endogen = "y", p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", trend = "restricted", r = 1,
                          error = "sv+covar", varsel = "bvs",
                          iterations = 10, burnin = 5)
  expect_error(
    add_priors(object, coef = list(v_i = 1, v_i_det = 1 / 10),
               coint = list(v_i = 1, p_tau_i = 1),
               sigma = list(mu = 0, v_i = 1 / 100, shape = 3, rate = 0.0001,
                            state_variance = 0.05, offset = 1e-8),
               varsel = list(inprior = 0.5, exclude_det = TRUE, covar = TRUE)),
    "one endogenous variable")
})
