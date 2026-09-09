export_to <- function(object, name) {
  folder <- file.path(tempdir(), name)
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  write_to_hdf5(object, folder = folder)
  folder
}

test_that("read_gvar_from_folder validates the folder", {
  skip_if_not_installed("hdf5r")

  expect_error(read_gvar_from_folder(file.path(tempdir(), "does-not-exist")),
               "does not exist")

  empty <- file.path(tempdir(), "bgvars-empty-export")
  unlink(empty, recursive = TRUE)
  dir.create(empty, recursive = TRUE)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)

  expect_error(read_gvar_from_folder(empty), "does not contain a model.h5")
})

test_that("a model without sub-models survives the round trip", {
  skip_if_not_installed("hdf5r")

  object <- gvar_object()
  folder <- export_to(object, "bgvars-read-plain")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)

  expect_s3_class(restored, "gvarmodel")
  expect_equal(names(restored), c("global", "weights", "submodels"))
  expect_null(restored[["submodels"]])
})

test_that("the global data survive the round trip", {
  skip_if_not_installed("hdf5r")

  object <- gvar_object()
  folder <- export_to(object, "bgvars-read-global")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)
  global <- restored[["global"]]

  expect_s3_class(global[["endogen"]], "ts")
  expect_equal(stats::tsp(global[["endogen"]]),
               stats::tsp(object[["global"]][["endogen"]]))
  expect_equal(dimnames(global[["endogen"]])[[2]],
               dimnames(object[["global"]][["endogen"]])[[2]])
  expect_equal(as.numeric(global[["endogen"]]),
               as.numeric(object[["global"]][["endogen"]]))

  expect_s3_class(global[["exogen"]], "ts")
  expect_equal(dimnames(global[["exogen"]])[[2]],
               dimnames(object[["global"]][["exogen"]])[[2]])

  expect_equal(global[["index"]], object[["global"]][["index"]])
})

test_that("a model without global data reads back without them", {
  skip_if_not_installed("hdf5r")

  submodel_data <- gvar_data()
  object <- create_gvarmodel(submodel_data = submodel_data)
  object <- add_weight_matrices(object, submodel_data, period = 3)

  folder <- export_to(object, "bgvars-read-noglobal")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)

  expect_null(restored[["global"]][["exogen"]])
  expect_s3_class(restored[["global"]][["endogen"]], "ts")
})

test_that("the weight matrices survive the round trip", {
  skip_if_not_installed("hdf5r")

  object <- gvar_object()
  folder <- export_to(object, "bgvars-read-weights")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)

  expect_equal(names(restored[["weights"]]), names(object[["weights"]]))
  for (i in names(object[["weights"]])) {
    expect_equal(restored[["weights"]][[i]], object[["weights"]][[i]],
                 ignore_attr = TRUE, info = i)
  }

  # A restored model can still be asked for a sub-model's weight matrix.
  expect_equal(get_weight_matrix(restored, "US"), get_weight_matrix(object, "US"))
})

test_that("sub-models come back named and in order", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10, p_endogen = 1:2)
  folder <- export_to(object, "bgvars-read-submodels")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)

  expect_equal(names(restored[["submodels"]]), names(object[["submodels"]]))
  for (i in names(object[["submodels"]])) {
    expect_s3_class(restored[["submodels"]][[i]], "modellist")
    expect_length(restored[["submodels"]][[i]], length(object[["submodels"]][[i]]))
    # The order the manifest lists them in is the order they were built in.
    expect_equal(vapply(restored[["submodels"]][[i]],
                        function(x) x[["model"]][["p_endogen"]], numeric(1)),
                 vapply(object[["submodels"]][[i]],
                        function(x) x[["model"]][["p_endogen"]], numeric(1)),
                 info = i)
  }
})

test_that("a sub-model keeps its specification, data and priors", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  folder <- export_to(object, "bgvars-read-submodel-content")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)

  original <- object[["submodels"]][["US"]][[1]]
  round_trip <- restored[["submodels"]][["US"]][[1]]

  expect_s3_class(round_trip, "bvarmodel")
  for (field in c("algorithm", "type", "k", "p_endogen", "k_exogen", "n",
                  "endogen", "exogen", "error", "iterations")) {
    expect_equal(round_trip[["model"]][[field]], original[["model"]][[field]],
                 info = field)
  }

  expect_equal(as.numeric(round_trip[["data"]][["train"]][["y"]]),
               as.numeric(original[["data"]][["train"]][["y"]]))
  expect_equal(round_trip[["data"]][["train"]][["z"]],
               original[["data"]][["train"]][["z"]], ignore_attr = TRUE)
  expect_equal(as.numeric(round_trip[["priors"]][["a"]][["mu"]]),
               as.numeric(original[["priors"]][["a"]][["mu"]]))
})

test_that("posterior draws written into the export are read back", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  folder <- export_to(object, "bgvars-read-posterior")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)
  draws <- restored[["submodels"]][["US"]][[1]][["posterior"]][["a"]][["coeffs"]]

  expect_false(is.null(draws))
  expect_equal(nrow(draws), 20)
  expect_equal(as.numeric(draws),
               as.numeric(object[["submodels"]][["US"]][[1]][["posterior"]][["a"]][["coeffs"]]))
})

test_that("a GVEC export comes back as a 'gvecmodel'", {
  skip_if_not_installed("hdf5r")

  object <- gvec_estimated(iterations = 20, burnin = 10)
  folder <- export_to(object, "bgvars-read-gvec")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)

  expect_s3_class(restored, "gvecmodel")
  expect_false(inherits(restored, "gvarmodel"))
  expect_s3_class(restored[["submodels"]][["US"]][[1]], "bvecmodel")
  expect_equal(restored[["submodels"]][["US"]][[1]][["model"]][["rank"]],
               object[["submodels"]][["US"]][[1]][["model"]][["rank"]])
})

test_that("sub-model identity comes from the manifest, not the directory tree", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  folder <- export_to(object, "bgvars-read-manifest")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  # Moving the whole export leaves every path different and every name the
  # same, which is the property reading the tree could not give.
  moved <- file.path(tempdir(), "bgvars-read-manifest-moved")
  unlink(moved, recursive = TRUE)
  file.rename(folder, moved)
  on.exit(unlink(moved, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(moved)
  expect_equal(names(restored[["submodels"]]), names(object[["submodels"]]))
})

test_that("a sub-model file the manifest names but cannot find is reported", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  folder <- export_to(object, "bgvars-read-missing")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  unlink(file.path(folder, "submodels", "JP", "001.h5"))

  expect_error(read_gvar_from_folder(folder), "submodels/JP/001.h5")
})

test_that("a restored model can be taken further", {
  skip_if_not_installed("hdf5r")

  object <- gvar_object()
  folder <- export_to(object, "bgvars-read-downstream")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  restored <- read_gvar_from_folder(folder)

  # Sub-models can be built on top of what came back, which is the test that
  # the parts a 'gvarmodel' is made of all arrived.
  restored <- add_submodels(restored,
                            endogen = c("y", "Dp"), p_endogen = 1,
                            exogen = c("y", "Dp"), p_exogen = 0,
                            global = "poil", s = 0,
                            iterations = 10, burnin = 10)

  expect_equal(names(restored[["submodels"]]), names(object[["weights"]]))
  expect_equal(restored[["submodels"]][["US"]][[1]][["data"]][["train"]][["y"]],
               create_varxsubmodel(object, submodel = "US",
                                   endogen = c("y", "Dp"), p_endogen = 1,
                                   exogen = c("y", "Dp"), p_exogen = 0,
                                   global = "poil", s = 0,
                                   iterations = 10, burnin = 10)[[1]][["data"]][["train"]][["y"]])
})
