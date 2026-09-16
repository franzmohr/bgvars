# Simulation and log-likelihoods of all models of a global model on several cores.

cores_gvec <- function() {
  object <- add_submodels(gvec_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", trend = "restricted", r = 0:1,
                          iterations = 20, burnin = 5)
  object <- add_priors(object, coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  add_seed(add_initial_values(object), 11)
}

cores_gvar <- function() {
  object <- add_submodels(gvar_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 1, deterministic = "const",
                          iterations = 20, burnin = 5)
  object <- add_priors(object, coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 1e-8))
  add_seed(add_initial_values(object), 12)
}

test_that("a gvecmodel draws the same on one core and on two", {
  skip_on_cran()
  object <- cores_gvec()
  one <- add_posterior_loglik(add_posterior_coefficients(object))
  two <- add_posterior_loglik(add_posterior_coefficients(object, cores = 2), cores = 2)

  expect_identical(names(two[["submodels"]]), names(one[["submodels"]]))
  expect_identical(lapply(two[["submodels"]], class), lapply(one[["submodels"]], class))
  expect_identical(lengths(two[["submodels"]]), lengths(one[["submodels"]]))
  for (s in names(one[["submodels"]])) {
    for (p in seq_along(one[["submodels"]][[s]])) {
      expect_equal(two[["submodels"]][[s]][[p]][["posterior"]],
                   one[["submodels"]][[s]][[p]][["posterior"]])
    }
  }
})

test_that("a gvarmodel draws the same on one core and on two", {
  skip_on_cran()
  object <- cores_gvar()
  one <- add_posterior_loglik(add_posterior_coefficients(object))
  two <- add_posterior_loglik(add_posterior_coefficients(object, cores = 2), cores = 2)

  for (s in names(one[["submodels"]])) {
    for (p in seq_along(one[["submodels"]][[s]])) {
      expect_equal(two[["submodels"]][[s]][[p]][["posterior"]],
                   one[["submodels"]][[s]][[p]][["posterior"]])
    }
  }
})

test_that("arguments for the samplers still reach every model", {
  object <- cores_gvec()
  marked <- add_posterior_coefficients(object, posterior_function = function(model) {
    model[["posterior"]] <- list(marker = TRUE)
    model
  })
  markers <- unlist(lapply(marked[["submodels"]], function(models) {
    vapply(models, function(model) isTRUE(model[["posterior"]][["marker"]]), logical(1))
  }))
  expect_true(all(markers))
  expect_length(markers, 6)
})
