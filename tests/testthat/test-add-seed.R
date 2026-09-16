# Seeds of the posterior simulation of the models of a global model.

seed_gvec <- function(r = 0:1) {
  object <- add_submodels(gvec_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", trend = "restricted", r = r,
                          iterations = 20, burnin = 5)
  add_priors(object,
             coef = list(v_i = 0),
             coint = list(v_i = 0, p_tau_i = 1),
             sigma = list(df = 3, scale = 0.0001))
}

seeds_of <- function(object) {
  unlist(lapply(object[["submodels"]], function(models) {
    vapply(models, function(model) model[["model"]][["seed"]], integer(1))
  }), use.names = FALSE)
}

test_that("add_initial_values() gives every model of a global model its own seed", {
  set.seed(1)
  object <- add_initial_values(seed_gvec())
  seeds <- seeds_of(object)

  # Three sub-models at ranks 0 and 1.
  expect_length(seeds, 6)
  expect_type(seeds, "integer")
  expect_identical(anyDuplicated(seeds), 0L)

  set.seed(1)
  expect_identical(seeds_of(add_initial_values(seed_gvec())), seeds)
})

test_that("add_seed() numbers the models through the sub-models in order", {
  object <- add_initial_values(seed_gvec())
  seeded <- add_seed(object, 50)

  expect_identical(seeds_of(seeded), 50L:55L)
  expect_identical(seeded[["submodels"]][[2]][[1]][["model"]][["seed"]], 52L)
  expect_identical(lapply(seeded[["submodels"]], class),
                   lapply(object[["submodels"]], class))
  expect_identical(names(seeded[["submodels"]]), names(object[["submodels"]]))
  expect_error(add_seed(object, -1), "whole number")
})

test_that("a seeded global model draws the same whatever the state of R's generator", {
  object <- add_seed(add_initial_values(seed_gvec(r = 1)), 7)

  set.seed(1)
  first <- add_posterior_coefficients(object)
  set.seed(2)
  again <- add_posterior_coefficients(object)

  for (submodel in names(object[["submodels"]])) {
    expect_identical(first[["submodels"]][[submodel]][[1]][["posterior"]],
                     again[["submodels"]][[submodel]][[1]][["posterior"]])
  }
})

test_that("add_seed() works for a gvarmodel", {
  object <- add_submodels(gvar_object(),
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 1,
                          deterministic = "const",
                          iterations = 20, burnin = 5)
  object <- add_priors(object, coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 1e-8))

  seeded <- add_seed(add_initial_values(object), 10)
  seeds <- seeds_of(seeded)
  expect_identical(seeds, seq(10L, length.out = length(seeds)))
})
