test_that("add_priors attaches priors on the cointegration space", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 10, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))

  expect_s3_class(object, "gvecmodel")
  for (i in names(object[["submodels"]])) {
    model <- object[["submodels"]][[i]][[1]]
    expect_true(all(c("a", "beta", "u_sigma") %in% names(model[["priors"]])))
  }
})

test_that("add_priors requires a specification for the cointegration space", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 10, burnin = 10)

  expect_error(add_priors(object, coef = list(v_i = 0),
                          sigma = list(df = 3, scale = 0.0001)))
})

test_that("add_initial_values attaches initial values to every sub-model", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 10, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)

  for (i in names(object[["submodels"]])) {
    initial <- object[["submodels"]][[i]][[1]][["initial"]]
    expect_true("beta" %in% names(initial))
    expect_true(all(is.finite(unlist(initial))))
  }
})

test_that("add_posterior_coefficients draws coefficients and beta", {
  iterations <- 40
  object <- gvec_estimated(iterations = iterations, burnin = 10)

  for (i in names(object[["submodels"]])) {
    model <- object[["submodels"]][[i]][[1]]
    posterior <- model[["posterior"]]
    expect_true(all(c("a", "beta", "u_sigma_inv") %in% names(posterior)))

    expect_equal(nrow(posterior[["a"]][["coeffs"]]), iterations)
    expect_true(all(is.finite(posterior[["a"]][["coeffs"]])))

    # One column per element of the cointegration matrix.
    expect_equal(nrow(posterior[["beta"]][["coeffs"]]), iterations)
    expect_equal(ncol(posterior[["beta"]][["coeffs"]]),
                 model[["model"]][["k_beta"]] * model[["model"]][["rank"]])
  }
})

test_that("a rank of zero drops the cointegration relation", {
  object <- gvec_estimated(iterations = 20, burnin = 10, r = 0)

  model <- object[["submodels"]][["US"]][[1]]
  expect_equal(model[["model"]][["rank"]], 0L)
  expect_true(all(is.finite(model[["posterior"]][["a"]][["coeffs"]])))
})

test_that("add_posterior_loglik adds one log-likelihood per draw and period", {
  iterations <- 20
  object <- gvec_estimated(iterations = iterations, burnin = 10)
  object <- add_posterior_loglik(object)

  for (i in names(object[["submodels"]])) {
    model <- object[["submodels"]][[i]][[1]]
    loglik <- model[["posterior"]][["loglik"]]
    expect_false(is.null(loglik))
    expect_equal(nrow(loglik), iterations)
    expect_equal(ncol(loglik), nrow(model[["data"]][["train"]][["y"]]))
  }
})
