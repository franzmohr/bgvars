# The thinning interval of the sampler of VEC and VAR sub-models.
#
# bvartools thins inside its samplers: a model that carries 'thin' runs
# burnin + iterations * thin draws and keeps the last of every 'thin', and the
# draws it keeps are labelled accordingly. A sub-model only has to carry the
# interval for its sampler to honour it, so what is tested here is that it does,
# and that the chain it produces is the unthinned chain with draws left out.

thin_submodel <- function(tvp, error, iterations, thin) {
  object <- gvec_object()
  models <- create_vecxsubmodel(object, submodel = "US",
                                endogen = c("y", "Dp"), p_endogen = 2,
                                exogen = c("y", "Dp"), p_exogen = 2,
                                const = "unrestricted", r = 1,
                                tvp = tvp, error = error,
                                iterations = iterations, burnin = 5, thin = thin)
  coef <- list(v_i = 1, v_i_det = 1 / 10)
  if (tvp) {
    coef <- c(coef, list(shape = 3, rate = 1e-8, rate_det = 1e-8))
  }
  sigma <- if (error == "wishart") {
    list(df = 3, scale = 0.0001)
  } else {
    list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
         state_variance = 0.05, offset = 1e-8)
  }
  coint <- if (tvp) list(rho = 0.999) else list(v_i = 0, p_tau_i = 1)
  add_priors(models, coef = coef, sigma = sigma, coint = coint)[[1]]
}

# Every block of draws in a posterior, named by where it sits.
draws_of <- function(posterior) {
  blocks <- list()
  for (name in names(posterior)) {
    element <- posterior[[name]]
    if (coda::is.mcmc(element)) {
      blocks[[name]] <- element
    } else if (is.list(element)) {
      for (inner in names(element)) {
        if (coda::is.mcmc(element[[inner]])) {
          blocks[[paste0(name, "$", inner)]] <- element[[inner]]
        }
      }
    }
  }
  blocks
}

test_that("create_vecxsubmodel carries the thinning interval of the sampler", {
  object <- gvec_object()
  build <- function(thin) {
    create_vecxsubmodel(object, submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 1,
                        exogen = c("y", "Dp"), p_exogen = 2, r = 0:1,
                        iterations = 10, burnin = 5, thin = thin)
  }

  for (model in build(1)) {
    expect_null(model[["model"]][["thin"]])
  }
  # Every model a call produces carries it, whatever its rank.
  for (model in build(4)) {
    expect_identical(model[["model"]][["thin"]], 4L)
  }

  for (thin in list(0, 2.5, c(2, 3), NA_real_, "2")) {
    expect_error(build(thin), "single positive integer")
  }
})

test_that("add_submodels passes the thinning interval to every sub-model", {
  object <- add_submodels(gvec_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2, r = 1,
                          iterations = 10, burnin = 5, thin = 3)

  for (submodel in object[["submodels"]]) {
    for (model in submodel) {
      expect_identical(model[["model"]][["thin"]], 3L)
    }
  }
})

test_that("a thinned chain is every thin-th draw of the unthinned one", {
  keep <- seq(3, 12, 3)

  for (spec in list(list(tvp = FALSE, error = "wishart"),
                    list(tvp = TRUE, error = "sv+covar"))) {
    set.seed(20260915)
    full <- add_posterior_coefficients(add_initial_values(
      thin_submodel(spec[["tvp"]], spec[["error"]], iterations = 12, thin = 1)))
    set.seed(20260915)
    thinned <- add_posterior_coefficients(add_initial_values(
      thin_submodel(spec[["tvp"]], spec[["error"]], iterations = 4, thin = 3)))

    before <- draws_of(full[["posterior"]])
    after <- draws_of(thinned[["posterior"]])
    expect_named(after, names(before))
    for (name in names(before)) {
      label <- paste(spec[["error"]], if (spec[["tvp"]]) "tvp", name)
      # Exactly: which draws are kept is the only thing thinning may change.
      expect_identical(c(as.matrix(after[[name]])),
                       c(as.matrix(before[[name]])[keep, , drop = FALSE]),
                       label = label)
      expect_equal(attr(after[[name]], "mcpar"), c(3, 12, 3), label = label)
    }
  }
})

# The same for the sub-models of a GVAR model.

thin_varx_submodel <- function(tvp, error, iterations, thin) {
  models <- create_varxsubmodel(gvar_object(), submodel = "US",
                                endogen = c("y", "Dp"), p_endogen = 1,
                                exogen = c("y", "Dp"), p_exogen = 1,
                                deterministic = "const",
                                tvp = tvp, error = error,
                                iterations = iterations, burnin = 5, thin = thin)
  coef <- list(v_i = 1, v_i_det = 1 / 10)
  if (tvp) {
    coef <- c(coef, list(shape = 3, rate = 1e-8, rate_det = 1e-8))
  }
  sigma <- if (error == "wishart") {
    list(df = 3, scale = 0.0001)
  } else {
    list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
         state_variance = 0.05, offset = 1e-8)
  }
  add_priors(models, coef = coef, sigma = sigma)[[1]]
}

test_that("create_varxsubmodel carries the thinning interval of the sampler", {
  object <- gvar_object()
  build <- function(thin) {
    create_varxsubmodel(object, submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 1:2,
                        exogen = c("y", "Dp"), p_exogen = 1,
                        iterations = 10, burnin = 5, thin = thin)
  }

  for (model in build(1)) {
    expect_null(model[["model"]][["thin"]])
  }
  for (model in build(4)) {
    expect_identical(model[["model"]][["thin"]], 4L)
  }

  for (thin in list(0, 2.5, c(2, 3), NA_real_, "2")) {
    expect_error(build(thin), "single positive integer")
  }
})

test_that("add_submodels passes the thinning interval to every GVAR sub-model", {
  object <- add_submodels(gvar_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          iterations = 10, burnin = 5, thin = 3)

  for (submodel in object[["submodels"]]) {
    for (model in submodel) {
      expect_identical(model[["model"]][["thin"]], 3L)
    }
  }
})

test_that("a thinned VARX chain is every thin-th draw of the unthinned one", {
  keep <- seq(3, 12, 3)

  for (spec in list(list(tvp = FALSE, error = "wishart"),
                    list(tvp = TRUE, error = "sv+covar"))) {
    set.seed(20260915)
    full <- add_posterior_coefficients(add_initial_values(
      thin_varx_submodel(spec[["tvp"]], spec[["error"]], iterations = 12, thin = 1)))
    set.seed(20260915)
    thinned <- add_posterior_coefficients(add_initial_values(
      thin_varx_submodel(spec[["tvp"]], spec[["error"]], iterations = 4, thin = 3)))

    before <- draws_of(full[["posterior"]])
    after <- draws_of(thinned[["posterior"]])
    expect_named(after, names(before))
    for (name in names(before)) {
      label <- paste(spec[["error"]], if (spec[["tvp"]]) "tvp", name)
      expect_identical(c(as.matrix(after[[name]])),
                       c(as.matrix(before[[name]])[keep, , drop = FALSE]),
                       label = label)
      expect_equal(attr(after[[name]], "mcpar"), c(3, 12, 3), label = label)
    }
  }
})
