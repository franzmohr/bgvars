# expected_model_size() for a global model and its export.
#
# Each model of a sub-model is sized by the method of bvartools, whose own tests
# hold the prediction against every sampler. What is checked here is the global
# side: that every model of every sub-model is counted once and labelled by
# where it is, that the prediction is what the estimation then stores, that an
# export is sized the same as the model it was written from, and that the
# warning of add_posterior_coefficients() is raised once for the whole global
# model rather than once per model.

size_gvar <- function() {
  add_priors(add_submodels(gvar_object(), endogen = c("y", "Dp"), p_endogen = 1:2,
                           exogen = c("y", "Dp"), p_exogen = 1,
                           deterministic = "const", iterations = 20, burnin = 5),
             coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8))
}

# The bytes of the coefficient draws of each model, predicted or stored.
predicted_bytes <- function(size) {
  size <- size[size[["step"]] == "add_posterior_coefficients", ]
  c(tapply(size[["bytes"]], size[["model"]], sum))
}

stored_bytes <- function(object) {
  models <- .submodel_models(object)[["models"]]
  names(models) <- .submodel_model_labels(rep(names(object[["submodels"]]),
                                              .submodel_models(object)[["sizes"]]))
  bytes <- vapply(models, function(model) {
    8 * sum(rapply(model[["posterior"]], length, how = "unlist"))
  }, numeric(1))
  bytes[order(names(bytes))]
}

test_that("every model of every sub-model is sized, and labelled by where it is", {
  object <- size_gvar()
  size <- expected_model_size(object)

  expect_s3_class(size, "modelsize")
  expect_setequal(unique(size[["model"]]),
                  c("global", "US/1", "US/2", "JP/1", "JP/2", "CA/1", "CA/2"))
  expect_equal(size[["bytes"]][size[["model"]] == "global"],
               as.numeric(utils::object.size(object[names(object) != "submodels"])))
  expect_output(print(size), "7 models")
})

test_that("the prediction is what the estimation stores", {
  object <- size_gvar()
  size <- expected_model_size(object)

  fitted <- add_posterior_coefficients(add_seed(add_initial_values(object), 42))
  expect_equal(predicted_bytes(size), stored_bytes(fitted))

  gvec <- add_submodels(gvec_object(), endogen = c("y", "Dp"), p_endogen = 1,
                        exogen = c("y", "Dp"), p_exogen = 1, const = "unrestricted",
                        r = 0:1, iterations = 20, burnin = 5)
  gvec <- add_priors(gvec, coef = list(v_i = 0), coint = list(v_i = 0, p_tau_i = 1),
                     sigma = list(df = 3, scale = 1e-4))
  size <- expected_model_size(gvec)
  fitted <- add_posterior_coefficients(add_seed(add_initial_values(gvec), 42))
  expect_equal(predicted_bytes(size), stored_bytes(fitted))
})

test_that("an export is sized as the model it was written from", {
  object <- size_gvar()
  folder <- file.path(tempdir(), "bgvars-export-size")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(add_initial_values(object), folder = folder, overwrite = TRUE, mc.cores = 1)
  export <- open_export(folder)

  in_memory <- expected_model_size(object)
  on_disk <- expected_model_size(export)
  expect_equal(predicted_bytes(on_disk), predicted_bytes(in_memory))
  expect_equal(on_disk[["bytes"]][on_disk[["model"]] == "global"],
               file.size(file.path(folder, "model.h5")))

  # The draws are not read, so an estimated export is sized the same.
  export <- add_posterior_coefficients(add_seed(export, 42))
  expect_equal(predicted_bytes(expected_model_size(export)), predicted_bytes(in_memory))
})

test_that("a global model too large to draw warns once, not once per model", {
  object <- add_initial_values(size_gvar())

  old <- options(bvartools.size_warning = 1000)
  on.exit(options(old), add = TRUE)

  warnings <- character(0)
  withCallingHandlers(add_posterior_coefficients(object),
                      warning = function(w) {
                        warnings <<- c(warnings, conditionMessage(w))
                        invokeRestart("muffleWarning")
                      })
  expect_length(warnings, 1)
  expect_match(warnings, "list of models")
})

test_that("a global model without sub-models is refused with what to do", {
  expect_error(expected_model_size(gvar_object()), "add_submodels")
})
