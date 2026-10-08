# time_variation_test() over a global model and over an export.
#
# The Bayes factors themselves are bvartools' and are tested there, against the
# ordinates the sampler stores. What is tested here is what only a global model
# has: that the results come back grouped per sub-model in the order of the
# sub-models, that an export gives one per stored model and is left untouched,
# and that the test reads rather than writes.

tvp_gvar_estimated <- function(iterations = 30, burnin = 10) {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          deterministic = "const", tvp = TRUE, error = "sv",
                          iterations = iterations, burnin = burnin)
  # Both blocks under the non-centred prior, so the test has both to report.
  object <- add_priors(object,
                       coef = list(v_i = 1, v_i_det = 1 / 10, omega_v = 0.001),
                       sigma = list(mu = 0, v_i = 1 / 100, omega_v = 0.1,
                                    state_variance = 0.05, offset = 1e-8))
  object <- add_initial_values(object)
  add_posterior_coefficients(object)
}

test_that("a global model is tested sub-model by sub-model", {

  object <- tvp_gvar_estimated()
  # A copy to compare against afterwards. Re-estimating would not do: the seed is
  # drawn from R's generator in add_initial_values(), so two estimations of the
  # same specification hold different draws.
  before <- object
  res <- time_variation_test(object)

  expect_type(res, "list")
  expect_named(res, names(object[["submodels"]]))
  for (submodel in names(res)) {
    expect_length(res[[submodel]], length(object[["submodels"]][[submodel]]))
    expect_s3_class(res[[submodel]][[1]], "bvartimevar")
    # The same numbers the single model gives.
    expect_equal(res[[submodel]][[1]],
                 time_variation_test(object[["submodels"]][[submodel]][[1]]))
  }

  # Both blocks are reported, which is what sigma$omega_v buys.
  expect_true(all(c("coefficients", "volatilities") %in% res[["US"]][[1]][["block"]]))

  # Arguments reach every model, and the test does not change the model.
  no_joint <- time_variation_test(object, joint = FALSE)
  expect_false(any(no_joint[["US"]][[1]][["term"]] == "(joint)"))
  expect_equal(object, before)
})

test_that("a global model without the non-centred prior is refused", {

  object <- gvar_estimated(iterations = 20, burnin = 5)
  expect_error(time_variation_test(object), "omega_v|non-centred")
})

test_that("an export is tested model by model and left as it was", {

  object <- tvp_gvar_estimated()
  folder <- file.path(tempdir(), "bgvars-timevar-export")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
  export <- read_export(folder)

  res <- time_variation_test(export)

  expect_type(res, "list")
  expect_length(res, nrow(export[["manifest"]]))
  expect_named(res, paste(export[["manifest"]][, "submodel"], "001", sep = "/"))
  expect_s3_class(res[[1]], "bvartimevar")

  # Reading leaves the manifest and the files exactly as they were.
  expect_equal(read_export(folder)[["manifest"]], export[["manifest"]])

  # One sub-model at a time.
  one <- time_variation_test(export, submodels = "US")
  expect_length(one, 1)
  expect_equal(one[[1]], res[["US/001"]])

  unlink(folder, recursive = TRUE)
})

test_that("map_submodels() reads without writing", {

  object <- gvar_estimated(iterations = 20, burnin = 5)
  folder <- file.path(tempdir(), "bgvars-map-read")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
  export <- read_export(folder)
  before <- export[["manifest"]]

  observations <- map_submodels(export, function(model) nrow(model$data$train$y),
                                write = FALSE)

  expect_type(observations, "list")
  expect_length(observations, nrow(before))
  expect_named(observations, paste(before[, "submodel"], "001", sep = "/"))
  expect_true(all(vapply(observations, is.numeric, logical(1))))
  expect_equal(read_export(folder)[["manifest"]], before)

  # A function that returns something other than a model is only allowed when
  # nothing is written.
  expect_error(map_submodels(export, function(model) nrow(model$data$train$y)),
               "rather than a model")
  expect_error(map_submodels(export, add_initial_values, write = NA),
               "must be TRUE or FALSE")

  unlink(folder, recursive = TRUE)
})
