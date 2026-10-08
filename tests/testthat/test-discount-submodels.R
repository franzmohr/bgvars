# Sub-models estimated by the two discounted models.
#
# What is checked is the specification the sub-model carries, the prior it gets,
# the manifest of an export of a grid of them, and the one thing they cannot do.
# The filter itself is BayesTS's and runs on the written files; bvartools checks
# the layout of those.

discount_gvec <- function(r = 0:2, ...) {
  object <- gvec_object(c("US", "CA"))
  add_submodels(object,
                endogen = c("y", "Dp"), p_endogen = 1,
                exogen = c("y", "Dp"), p_exogen = 2,
                const = "unrestricted", r = r,
                algorithm = "discount",
                iterations = 20, burnin = 0, thin = 1, ...)
}


test_that("a sub-model can be a discounted model", {

  object <- discount_gvec(r = 1, delta_beta = 0.98, delta_sigma = 0.95)
  specification <- object[["submodels"]][["US"]][[1]][["model"]]

  expect_equal(specification[["algorithm"]], "VecTvpDiscount")
  expect_equal(specification[["delta_beta"]], 0.98)
  expect_equal(specification[["delta_sigma"]], 0.95)
  expect_true(specification[["tvp"]])
  expect_true(bvartools::is_discount_model(object[["submodels"]][["US"]][[1]]))

  expect_error(discount_gvec(r = 1, delta_beta = 1.5), "delta_beta")
  expect_error(add_submodels(gvec_object(c("US", "CA")),
                             endogen = c("y", "Dp"), p_endogen = 1,
                             exogen = c("y", "Dp"), p_exogen = 2,
                             const = "unrestricted", r = 1,
                             algorithm = "discount", error = "sv",
                             iterations = 20, burnin = 0, thin = 1),
               "wishart")
})


test_that("a grid over ranks and discounts gives one sub-model per combination", {

  object <- discount_gvec(r = 0:2, delta_beta = c(0.98, 1))

  for (submodel in names(object[["submodels"]])) {
    models <- object[["submodels"]][[submodel]]
    expect_length(models, 6)
    grid <- t(vapply(models, function(m) {
      c(m[["model"]][["rank"]], m[["model"]][["delta_beta"]])
    }, numeric(2)))
    expect_equal(nrow(unique(grid)), 6)
  }
})


test_that("the prior of a discounted sub-model is the matrix normal one", {

  object <- discount_gvec(r = 1)
  object <- add_priors(object, coef = list(v_i = 1, v_i_det = 1 / 10),
                       sigma = list(df = "k", scale = 1))

  model <- object[["submodels"]][["US"]][[1]]
  k <- model[["model"]][["k"]]
  n_design <- model[["model"]][["rank"]] + ncol(model[["data"]][["train"]][["x"]])

  prior <- model[["priors"]][["a"]]
  expect_equal(dim(prior[["mean"]]), c(n_design, k))
  expect_equal(dim(prior[["cov"]]), c(n_design, n_design))
  expect_null(prior[["v_inv"]])
  # No cointegration space prior: the space is conditioned on rather than drawn,
  # and 'coint' is never given.
  expect_null(model[["priors"]][["beta"]])

  object <- add_initial_values(object)
  model <- object[["submodels"]][["US"]][[1]]
  expect_equal(length(model[["initial"]][["beta"]]),
               ncol(model[["data"]][["train"]][["w"]]) * 1L)
})


test_that("an export of discounted sub-models records what they differ in", {

  object <- discount_gvec(r = 0:2, delta_beta = c(0.98, 1))
  object <- add_priors(object, coef = list(v_i = 1, v_i_det = 1 / 10),
                       sigma = list(df = "k", scale = 1))
  object <- add_initial_values(object)

  folder <- file.path(tempdir(), "bgvars-discount-export")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  export <- read_export(folder)
  manifest <- export[["manifest"]]
  expect_equal(nrow(manifest), 12L)
  expect_true(all(c("delta_beta", "delta_sigma") %in% names(manifest)))
  expect_setequal(manifest[["delta_beta"]], c(0.98, 1))
  # One is what a model without discounts is, rather than missing.
  expect_true(all(manifest[["delta_sigma"]] == 1))

  # The SUR matrix is not written for these, so a file is smaller than the
  # sampler's; what matters here is that the export is still readable.
  expect_length(submodel_directories(export), 2L)
})


test_that("a closed form cannot be stacked into a global model", {

  object <- discount_gvec(r = 1)
  object <- add_priors(object, coef = list(v_i = 1, v_i_det = 1 / 10),
                       sigma = list(df = "k", scale = 1))
  object <- add_initial_values(object)

  # A global model is built by stacking one draw of every sub-model at a time,
  # and a closed form has none: joining one draw per period would look like a
  # sampled path and is not one.
  expect_error(submodels_to_gvar(object), "closed-form posterior")

  folder <- file.path(tempdir(), "bgvars-discount-solve")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  # From the manifest, before a file is read.
  expect_error(submodels_to_gvar(read_export(folder)), "closed-form posterior")
})


test_that("a discounted sub-model carries no SUR matrix", {

  object <- discount_gvec(r = 0:1)
  for (model in object[["submodels"]][["US"]]) {
    expect_null(model[["data"]][["train"]][["z"]])
    expect_false(is.null(model[["data"]][["train"]][["x"]]))
  }
})
