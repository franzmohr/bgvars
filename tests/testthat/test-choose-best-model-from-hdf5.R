# Choosing one model per sub-model from an export, and reading the choice back.

export_grid <- function(name = "bgvars-selection") {
  object <- gvar_estimated_grid(iterations = 30, burnin = 10)
  folder <- file.path(tempdir(), name)
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  write_to_hdf5(object, folder = folder)
  list(object = object, folder = folder)
}

test_that("choose_best_model_from_hdf5 validates the folder", {
  skip_if_not_installed("hdf5r")

  expect_error(
    choose_best_model_from_hdf5(file.path(tempdir(), "does-not-exist"), "BIC"),
    "does not exist")

  empty <- file.path(tempdir(), "bgvars-selection-empty")
  unlink(empty, recursive = TRUE)
  dir.create(empty, recursive = TRUE)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)

  expect_error(choose_best_model_from_hdf5(empty, "BIC"),
               "does not contain a model.h5")
})

test_that("choose_best_model_from_hdf5 returns one position per sub-model", {
  skip_if_not_installed("hdf5r")

  export <- export_grid()
  on.exit(unlink(export[["folder"]], recursive = TRUE), add = TRUE)

  result <- choose_best_model_from_hdf5(export[["folder"]], criterion = "BIC")

  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("submodel", "position"))
  expect_equal(result[, "submodel"], names(export[["object"]][["submodels"]]))
  expect_true(all(result[, "position"] %in% 1:2))
})

test_that("reading a whole sub-model gives the same ranking", {
  skip_if_not_installed("hdf5r")

  export <- export_grid("bgvars-selection-agree")
  on.exit(unlink(export[["folder"]], recursive = TRUE), add = TRUE)

  from_folder <- choose_best_model_from_hdf5(export[["folder"]],
                                             criterion = "BIC")

  # The route the vignette shows for a single sub-model: read every model of it
  # and rank them with bvartools. Reading only the inputs of the criteria has
  # to give the same answer as reading the models in full.
  for (submodel in from_folder[, "submodel"]) {
    models <- read_submodel_from_folder(export[["folder"]], submodel)
    expect_equal(choose_best_model(selection_criteria(models),
                                   criterion = "BIC"),
                 from_folder[from_folder[, "submodel"] == submodel, "position"],
                 info = submodel)
  }
})

test_that("read_submodel_from_folder validates its arguments", {
  skip_if_not_installed("hdf5r")

  export <- export_grid("bgvars-selection-read")
  on.exit(unlink(export[["folder"]], recursive = TRUE), add = TRUE)

  expect_error(read_submodel_from_folder(file.path(tempdir(), "nope"), "US"),
               "does not exist")
  expect_error(read_submodel_from_folder(export[["folder"]], c("US", "JP")),
               "must have length 1")
  expect_error(read_submodel_from_folder(export[["folder"]], 1),
               "must be of class 'character'")
  expect_error(read_submodel_from_folder(export[["folder"]], "XX"),
               "is not available")
})

test_that("read_gvar_from_folder keeps the models a selection names", {
  skip_if_not_installed("hdf5r")

  export <- export_grid("bgvars-selection-restore")
  on.exit(unlink(export[["folder"]], recursive = TRUE), add = TRUE)

  best <- choose_best_model_from_hdf5(export[["folder"]], criterion = "BIC")
  restored <- read_gvar_from_folder(export[["folder"]], submodels = best)

  expect_equal(vapply(restored[["submodels"]], length, numeric(1)),
               stats::setNames(rep(1, nrow(best)), best[, "submodel"]))

  # The position of a model is its position among the models of its sub-model,
  # so the lag order that comes back is the one that was selected.
  expected <- vapply(seq_len(nrow(best)), function(i) {
    submodel <- best[i, "submodel"]
    models <- export[["object"]][["submodels"]][[submodel]]
    models[[best[i, "position"]]][["model"]][["p_endogen"]]
  }, numeric(1))

  expect_equal(unname(vapply(restored[["submodels"]],
                             function(x) x[[1]][["model"]][["p_endogen"]],
                             numeric(1))),
               expected)
})

test_that("a selected global model can be solved", {
  skip_if_not_installed("hdf5r")

  export <- export_grid("bgvars-selection-solve")
  on.exit(unlink(export[["folder"]], recursive = TRUE), add = TRUE)

  best <- choose_best_model_from_hdf5(export[["folder"]], criterion = "BIC")
  restored <- read_gvar_from_folder(export[["folder"]], submodels = best)

  gvar <- submodels_to_gvar(restored)

  expect_s3_class(gvar, "bvarmodel")
  expect_equal(gvar[["model"]][["k"]], 6L)
})

test_that("a selection that does not fit the export is reported", {
  skip_if_not_installed("hdf5r")

  export <- export_grid("bgvars-selection-invalid")
  on.exit(unlink(export[["folder"]], recursive = TRUE), add = TRUE)

  folder <- export[["folder"]]

  expect_error(
    read_gvar_from_folder(folder,
                          submodels = data.frame(submodel = "XX", position = 1)),
    "does not contain the sub-model")

  expect_error(
    read_gvar_from_folder(folder,
                          submodels = data.frame(submodel = "US", position = 99)),
    "position 99 does not exist")

  expect_error(
    read_gvar_from_folder(folder, submodels = data.frame(submodel = "US")),
    "must have the column")
})
