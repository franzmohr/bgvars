test_that("write_to_hdf5 requires an existing folder", {
  skip_if_not_installed("hdf5r")

  object <- gvar_object()
  folder <- file.path(tempdir(), "does-not-exist")

  expect_error(write_to_hdf5(object, folder = folder), "does not exist")
})

test_that("write_to_hdf5 exports a GVAR model to a folder", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)

  folder <- file.path(tempdir(), "bgvars-gvar-export")
  unlink(folder, recursive = TRUE)
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  write_to_hdf5(object, folder = folder)

  expect_true(file.exists(file.path(folder, "global.h5")))
  expect_true(file.exists(file.path(folder, "weights.h5")))
  expect_setequal(list.files(file.path(folder, "submodels")),
                  names(object[["submodels"]]))

  global <- hdf5r::h5file(file.path(folder, "global.h5"), mode = "r")
  on.exit(global$close_all(), add = TRUE)
  expect_setequal(names(global), c("endogen", "exogen", "index"))
  expect_equal(hdf5r::h5attr(global[["endogen"]], "variables"),
               dimnames(object[["global"]][["endogen"]])[[2]])
  expect_equal(hdf5r::h5attr(global[["endogen"]], "tsp"),
               as.numeric(stats::tsp(object[["global"]][["endogen"]])))
})

test_that("write_to_hdf5 exports a GVEC model to a folder", {
  skip_if_not_installed("hdf5r")

  object <- gvec_estimated(iterations = 20, burnin = 10)

  folder <- file.path(tempdir(), "bgvars-gvec-export")
  unlink(folder, recursive = TRUE)
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  write_to_hdf5(object, folder = folder)

  expect_true(file.exists(file.path(folder, "global.h5")))
  expect_true(file.exists(file.path(folder, "weights.h5")))
  expect_setequal(list.files(file.path(folder, "submodels")),
                  names(object[["submodels"]]))
})
