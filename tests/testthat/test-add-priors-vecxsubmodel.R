# Priors of a VECX sub-model that are specified by name rather than by value.

vecx_submodel <- function() {
  object <- gvec_object()
  models <- suppressMessages(
    create_vecxsubmodel(object, submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 2,
                        exogen = c("y", "Dp"), p_exogen = 2,
                        const = "unrestricted", r = 1,
                        iterations = 10, burnin = 10))
  models[[1]]
}

test_that("sigma$df can be given in terms of k", {
  model <- vecx_submodel()

  result <- add_priors(model,
                       coef = list(v_i = 1),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = "k + 3", scale = 0.0001))

  # k = 2 endogenous variables plus 3. The rank is not added: BayesTS adds it
  # to the posterior degrees of freedom itself.
  expect_equal(result[["priors"]][["u_sigma"]][["df"]], 2 + 3)
})

test_that("sigma$df and sigma$shape are stored as given for any rank", {
  submodel <- function(r, error) {
    suppressMessages(
      create_vecxsubmodel(gvec_object(), submodel = "US",
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", r = r, error = error,
                          iterations = 10, burnin = 10))[[1]]
  }

  for (r in 1:2) {
    result <- add_priors(submodel(r, "wishart"), coef = list(v_i = 1),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = 3, scale = 0.0001))
    expect_equal(result[["priors"]][["u_sigma"]][["type"]], "wishart")
    expect_equal(result[["priors"]][["u_sigma"]][["df"]], 3)

    result <- add_priors(submodel(r, "gamma"), coef = list(v_i = 1),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(shape = 3, rate = 0.0001))
    expect_equal(result[["priors"]][["u_sigma"]][["type"]], "gamma")
    expect_equal(c(result[["priors"]][["u_sigma"]][["shape"]]), c(3, 3))
  }
})

test_that("coef$const sets the prior mean of the unrestricted constant", {
  model <- vecx_submodel()
  k <- model[["model"]][["k"]]
  r <- model[["model"]][["rank"]]
  column <- r + which(dimnames(model[["data"]][["train"]][["x"]])[[2]] == "const")

  prior_mean <- function(const) {
    result <- add_priors(model,
                         coef = list(v_i = 1, const = const),
                         coint = list(v_i = 0, p_tau_i = 1),
                         sigma = list(df = 3, scale = 0.0001))
    matrix(result[["priors"]][["a"]][["mu"]], k)
  }

  mu <- prior_mean(0.5)
  expect_equal(mu[, column], rep(0.5, k))
  expect_true(all(mu[, -column] == 0))

  expect_equal(prior_mean(c(0.1, 0.2))[, column], c(0.1, 0.2))

  expect_equal(prior_mean("mean")[, column],
               unname(colMeans(model[["data"]][["train"]][["y"]])))
  expect_equal(prior_mean("first")[, column],
               unname(model[["data"]][["train"]][["y"]][1, ]))

  expect_error(prior_mean(c(0.1, 0.2, 0.3)), "coef[$]const")
})
