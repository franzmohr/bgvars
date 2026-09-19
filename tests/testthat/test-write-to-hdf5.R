export_folder <- function(name) {
  folder <- file.path(tempdir(), name)
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  folder
}

test_that("write_to_hdf5 requires an existing folder", {
  skip_if_not_installed("hdf5r")

  object <- gvar_object()
  expect_error(write_to_hdf5(object, folder = file.path(tempdir(), "does-not-exist")),
               "does not exist")
})

test_that("write_to_hdf5 writes one model file and one file per sub-model", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10, p_endogen = 1:2)

  folder <- export_folder("bgvars-gvar-export")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  expect_equal(write_to_hdf5(object, folder = folder), folder)

  expect_true(file.exists(file.path(folder, "model.h5")))
  # The global data and the weights share one file now.
  expect_false(file.exists(file.path(folder, "global.h5")))
  expect_false(file.exists(file.path(folder, "weights.h5")))

  expect_setequal(list.files(file.path(folder, "submodels")),
                  names(object[["submodels"]]))
  # Two lag orders were requested, so every sub-model contributes two models,
  # named by position rather than by specification.
  for (i in names(object[["submodels"]])) {
    expect_equal(list.files(file.path(folder, "submodels", i)),
                 c("001.h5", "002.h5"))
  }
})

test_that("model.h5 holds the global data, the weights and the manifest", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)

  folder <- export_folder("bgvars-model-file")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  h5 <- hdf5r::h5file(file.path(folder, "model.h5"), mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  expect_setequal(names(h5), c("global", "weights", "submodels"))
  expect_equal(hdf5r::h5attr(h5, "rclass"), class(object))

  expect_setequal(names(h5[["global"]]),
                  c("endogen", "exogen", "deterministic", "index"))
  expect_equal(hdf5r::h5attr(h5[["global"]][["endogen"]], "variables"),
               dimnames(object[["global"]][["endogen"]])[[2]])
  expect_equal(hdf5r::h5attr(h5[["global"]][["endogen"]], "tsp"),
               as.numeric(stats::tsp(object[["global"]][["endogen"]])))

  expect_setequal(names(h5[["weights"]]), names(object[["weights"]]))
  expect_equal(h5[["weights"]][["US"]]$dims, dim(object[["weights"]][["US"]]))
  expect_equal(h5[["weights"]][["US"]]$read(), object[["weights"]][["US"]])
})

test_that("the manifest describes every sub-model file that was written", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10, p_endogen = 1:2)

  folder <- export_folder("bgvars-manifest")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  h5 <- hdf5r::h5file(file.path(folder, "model.h5"), mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  manifest <- h5[["submodels"]]$read()

  expect_s3_class(manifest, "data.frame")
  expect_equal(nrow(manifest),
               sum(vapply(object[["submodels"]], length, numeric(1))))
  expect_setequal(manifest[, "submodel"], names(object[["submodels"]]))

  # Every row points at a file that is really there ...
  for (i in manifest[, "file"]) {
    expect_true(file.exists(file.path(folder, i)))
  }
  # ... and every file that is there has a row.
  written <- list.files(file.path(folder, "submodels"), recursive = TRUE)
  expect_setequal(paste0("submodels/", written), manifest[, "file"])

  # The specification the file names used to carry.
  expect_true(all(manifest[, "algorithm"] == "VarNormalWishart"))
  expect_true(all(manifest[, "type"] == "VARX"))
  expect_setequal(manifest[manifest[, "submodel"] == "US", "p_endogen"], 1:2)
  expect_true(all(manifest[, "iterations"] == 20))
  # A model that keeps every draw is recorded with a thinning interval of one.
  expect_true(all(manifest[, "thin"] == 1L))
  # A VARX has no cointegration rank.
  expect_true(all(is.na(manifest[, "rank"])))
})

test_that("the manifest records the cointegration rank of a GVEC model", {
  skip_if_not_installed("hdf5r")

  object <- gvec_estimated(iterations = 20, burnin = 10, r = 0:1)

  folder <- export_folder("bgvars-gvec-export")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  h5 <- hdf5r::h5file(file.path(folder, "model.h5"), mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  manifest <- h5[["submodels"]]$read()

  expect_true(all(manifest[, "type"] == "VECX"))
  expect_setequal(manifest[manifest[, "submodel"] == "US", "rank"], 0:1)
  expect_equal(hdf5r::h5attr(h5, "rclass"), class(object))
})

test_that("the manifest records the thinning interval of the sampler", {
  skip_if_not_installed("hdf5r")

  object <- gvec_estimated(iterations = 10, burnin = 5, thin = 3)

  folder <- export_folder("bgvars-thin-export")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  h5 <- hdf5r::h5file(file.path(folder, "model.h5"), mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  manifest <- h5[["submodels"]]$read()

  expect_true(all(manifest[, "thin"] == 3L))
  # 'iterations' stays the number of draws kept.
  expect_true(all(manifest[, "iterations"] == 10L))
})

test_that("a sub-model file carries the model, its data and its priors", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)

  folder <- export_folder("bgvars-submodel-file")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  h5 <- hdf5r::h5file(file.path(folder, "submodels", "US", "001.h5"), mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  expect_true(all(c("model", "data", "priors") %in% names(h5)))
  # The algorithm attribute is what an external sampler dispatches on.
  expect_equal(hdf5r::h5attr(h5[["model"]], "algorithm"), "VarNormalWishart")
  expect_equal(h5[["data"]][["train"]][["y"]]$dims,
               dim(object[["submodels"]][["US"]][[1]][["data"]][["train"]][["y"]]))
})

test_that("write_to_hdf5 refuses to overwrite an export unless asked", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)

  folder <- export_folder("bgvars-overwrite")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  expect_error(write_to_hdf5(object, folder = folder), "already exists")
  expect_no_error(write_to_hdf5(object, folder = folder, overwrite = TRUE))

  # Replaced, not added to.
  expect_equal(list.files(file.path(folder, "submodels", "US")), "001.h5")
})

test_that("a failure part way through the export is reported", {
  skip_if_not_installed("hdf5r")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  # Something bvartools has no method for, so that the export fails after
  # model.h5 has been written. The function used to wrap its whole body in
  # try(), which turned this into a silent return and a folder that looked
  # finished.
  object[["submodels"]][["JP"]][[1]] <- list(model = list(algorithm = "nonsense"))

  folder <- export_folder("bgvars-failure")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  expect_error(write_to_hdf5(object, folder = folder))
})

test_that("weight matrices are stored compressed", {
  skip_if_not_installed("hdf5r")

  object <- gvar_object()

  folder <- export_folder("bgvars-compression")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  h5 <- hdf5r::h5file(file.path(folder, "model.h5"), mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  weights <- h5[["weights"]][["US"]]
  uncompressed <- prod(dim(object[["weights"]][["US"]])) * 8
  expect_lt(weights$get_storage_size(), uncompressed)
  # The matrices are mostly zeros, so the saving should be substantial.
  expect_lt(weights$get_storage_size(), uncompressed / 4)
})

test_that(".hdf5_chunk_dims spans the full width of a matrix", {
  x <- matrix(0, 2000, 175)

  chunk <- bgvars:::.hdf5_chunk_dims(x)
  expect_equal(chunk[2], ncol(x))
  expect_lte(chunk[1], nrow(x))

  # A matrix small enough to fit the budget is one chunk.
  expect_equal(bgvars:::.hdf5_chunk_dims(matrix(0, 10, 10)), c(10, 10))

  # A wide matrix is split by rows rather than given an oversized chunk.
  wide <- matrix(0, 5000, 5000)
  chunk_wide <- bgvars:::.hdf5_chunk_dims(wide, max_bytes = 8e6)
  expect_equal(chunk_wide[2], 5000)
  expect_lt(chunk_wide[1], 5000)
  expect_lte(prod(chunk_wide) * 8, 8e6)
})

test_that("an export leaves nothing open in the files it wrote", {
  skip_if_not_installed("hdf5r")

  # HDF5 keeps a file open for as long as anything in it is, even once the
  # handle of the file is closed, and on Windows the file is then locked
  # against every other process -- BayesTS among them. Counting straight after
  # the export, through a second handle that is not counted itself, makes this
  # independent of when the garbage collector would release what was left. It
  # guards the export as a whole: today its close_all() would release a handle
  # left open, and this fails if the files come to be closed some other way.
  open_in_h5_file <- function(path) {
    h5 <- hdf5r::H5File$new(path, mode = "r")
    on.exit(h5$close(), add = TRUE)
    h5$get_obj_count() - 1
  }

  object <- gvar_estimated(iterations = 20, burnin = 10)

  folder <- export_folder("bgvars-open-handles")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  expect_equal(open_in_h5_file(file.path(folder, "model.h5")), 0)
  for (i in names(object[["submodels"]])) {
    path <- file.path(folder, "submodels", i, "001.h5")
    expect_equal(open_in_h5_file(path), 0, info = i)
  }
})

# --- Appending one batch of sub-models at a time -------------------------------
#
# A global model is written batch by batch, which is what keeps a run that
# cannot hold every sub-model's draws from having to.

append_object <- function() {
  object <- gvar_object()
  object <- add_submodels(object, endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          deterministic = "const", iterations = 20, burnin = 5)
  object <- add_priors(object, coef = list(v_i = 1),
                       sigma = list(df = 3, scale = 1e-8))
  add_seed(add_initial_values(object), 31)
}

batch_of <- function(object, submodels) {
  object[["submodels"]] <- object[["submodels"]][submodels]
  object
}

test_that("sub-models written in batches give the export a single write gives", {
  object <- add_posterior_coefficients(append_object())
  submodels <- names(object[["submodels"]])
  expect_gte(length(submodels), 3)

  whole <- file.path(tempdir(), "bgvars-append-whole")
  batched <- file.path(tempdir(), "bgvars-append-batched")
  for (folder in c(whole, batched)) {
    unlink(folder, recursive = TRUE)
    dir.create(folder, recursive = TRUE)
  }

  write_to_hdf5(object, folder = whole, overwrite = TRUE, mc.cores = 1)
  first <- submodels[1]
  rest <- submodels[-1]
  write_to_hdf5(batch_of(object, first), folder = batched, overwrite = TRUE, mc.cores = 1)
  write_to_hdf5(batch_of(object, rest), folder = batched, append = TRUE, mc.cores = 1)

  manifest_of <- function(folder) {
    h5 <- hdf5r::h5file(file.path(folder, "model.h5"), mode = "r")
    on.exit(h5$close_all(), add = TRUE)
    as.data.frame(h5[["submodels"]]$read(), stringsAsFactors = FALSE)
  }
  one <- manifest_of(whole)
  many <- manifest_of(batched)
  expect_setequal(many[["submodel"]], one[["submodel"]])
  expect_identical(nrow(many), nrow(one))
  expect_equal(many[order(many[["submodel"]]), ], one[order(one[["submodel"]]), ],
               ignore_attr = TRUE)

  # The draws of every sub-model survive the batching.
  back <- read_gvar_from_folder(batched)
  for (submodel in submodels) {
    expect_equal(back[["submodels"]][[submodel]][[1]][["posterior"]][["a"]][["coeffs"]],
                 object[["submodels"]][[submodel]][[1]][["posterior"]][["a"]][["coeffs"]])
  }
  # The global data and the weights come from the first batch and describe the
  # whole model.
  expect_equal(back[["global"]][["index"]], object[["global"]][["index"]])
  expect_equal(names(back[["weights"]]), names(object[["weights"]]))
})

test_that("appending a sub-model that is there replaces it", {
  object <- add_posterior_coefficients(append_object())
  submodels <- names(object[["submodels"]])
  folder <- file.path(tempdir(), "bgvars-append-replace")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)

  write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)

  changed <- batch_of(object, submodels[1])
  coeffs <- changed[["submodels"]][[1]][[1]][["posterior"]][["a"]][["coeffs"]]
  changed[["submodels"]][[1]][[1]][["posterior"]][["a"]][["coeffs"]][] <- coeffs * 0 + 1
  write_to_hdf5(changed, folder = folder, append = TRUE, mc.cores = 1)

  back <- read_gvar_from_folder(folder)
  expect_length(back[["submodels"]], length(submodels))
  expect_true(all(unclass(back[["submodels"]][[submodels[1]]][[1]][["posterior"]][["a"]][["coeffs"]]) == 1))
  expect_equal(back[["submodels"]][[submodels[2]]][[1]][["posterior"]][["a"]][["coeffs"]],
               object[["submodels"]][[submodels[2]]][[1]][["posterior"]][["a"]][["coeffs"]])
})

test_that("appending to a folder without an export is refused", {
  object <- add_posterior_coefficients(append_object())
  folder <- file.path(tempdir(), "bgvars-append-empty")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  expect_error(write_to_hdf5(object, folder = folder, append = TRUE, mc.cores = 1),
               "no export to append to")
})
