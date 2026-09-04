test_that("add_priors attaches priors to every sub-model", {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 0,
                          global = "poil", s = 0,
                          iterations = 10, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 0.0001))

  expect_s3_class(object, "gvarmodel")
  for (i in names(object[["submodels"]])) {
    model <- object[["submodels"]][[i]][[1]]
    expect_true("priors" %in% names(model))
    expect_true(all(c("a", "u_sigma") %in% names(model[["priors"]])))

    # One prior mean and one prior precision per coefficient.
    n_coef <- ncol(model[["data"]][["train"]][["z"]])
    expect_equal(model[["priors"]][["a"]][["type"]], "normal")
    expect_length(model[["priors"]][["a"]][["mu"]], n_coef)
    expect_equal(dim(model[["priors"]][["a"]][["v_inv"]]), c(n_coef, n_coef))
    expect_true(all(diag(model[["priors"]][["a"]][["v_inv"]]) == 1))
    expect_equal(model[["priors"]][["u_sigma"]][["type"]], "wishart")
    expect_equal(model[["priors"]][["u_sigma"]][["df"]], 3L)
  }
})

test_that("add_initial_values attaches initial values to every sub-model", {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 0,
                          global = "poil", s = 0,
                          iterations = 10, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)

  for (i in names(object[["submodels"]])) {
    model <- object[["submodels"]][[i]][[1]]
    initial <- model[["initial"]]
    expect_true(all(c("a", "u_sigma_inv") %in% names(initial)))
    expect_length(initial[["a"]], ncol(model[["data"]][["train"]][["z"]]))
    expect_true(all(is.finite(initial[["a"]])))
    k <- model[["model"]][["k"]]
    expect_equal(dim(initial[["u_sigma_inv"]]), c(k, k))
    # A precision matrix must be symmetric.
    expect_equal(initial[["u_sigma_inv"]], t(initial[["u_sigma_inv"]]))
  }
})

test_that("add_posterior_coefficients draws the requested number of samples", {
  iterations <- 40
  object <- gvar_estimated(iterations = iterations, burnin = 10)

  for (i in names(object[["submodels"]])) {
    model <- object[["submodels"]][[i]][[1]]
    posterior <- model[["posterior"]]
    expect_true(all(c("a", "u_sigma_inv") %in% names(posterior)))

    draws_a <- posterior[["a"]][["coeffs"]]
    expect_s3_class(draws_a, "mcmc")
    expect_equal(nrow(draws_a), iterations)
    expect_equal(ncol(draws_a), ncol(model[["data"]][["train"]][["z"]]))
    expect_true(all(is.finite(draws_a)))

    draws_sigma <- posterior[["u_sigma_inv"]][["coeffs"]]
    expect_equal(nrow(draws_sigma), iterations)
    expect_equal(ncol(draws_sigma), model[["model"]][["k"]]^2)
  }
})

test_that("posterior simulation leaves the model set-up untouched", {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 0,
                          global = "poil", s = 0,
                          iterations = 20, burnin = 10)
  before <- object[["submodels"]][["US"]][[1]]

  object <- add_priors(object, coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)
  after <- object[["submodels"]][["US"]][[1]]

  expect_equal(after[["model"]], before[["model"]])
  expect_equal(after[["data"]], before[["data"]])
})

test_that("add_posterior_loglik adds one log-likelihood per draw and period", {
  iterations <- 20
  object <- gvar_estimated(iterations = iterations, burnin = 10)
  object <- add_posterior_loglik(object)

  for (i in names(object[["submodels"]])) {
    model <- object[["submodels"]][[i]][[1]]
    loglik <- model[["posterior"]][["loglik"]]
    expect_false(is.null(loglik))
    expect_equal(nrow(loglik), iterations)
    expect_equal(ncol(loglik), nrow(model[["data"]][["train"]][["y"]]))
    expect_true(all(is.finite(loglik)))
  }
})

test_that("the estimated model can be summarised by selection criteria", {
  object <- gvar_estimated(iterations = 20, burnin = 10, p_endogen = 1:2)
  object <- add_posterior_loglik(object)

  criteria <- bvartools::selection_criteria(object[["submodels"]][["US"]])

  expect_length(criteria, 2)
  for (i in criteria) {
    expect_true(all(c("LL", "AIC", "BIC", "HQ") %in% names(i)))
  }
})
