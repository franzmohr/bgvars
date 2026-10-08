# The persistence profiles of a solved global model. What is checked here is the
# part only a global model has: that the cointegrating vectors of the sub-models
# arrive in global coordinates and that the relations can be followed whole. The
# expression itself belongs to bvartools' test-persistence_profiles.R.

test_that("a solved GVEC carries its cointegrating vectors in global coordinates", {
  object <- gvec_estimated(iterations = 30, burnin = 10, r = 1)
  gvar <- suppressMessages(submodels_to_gvar(object))

  # One relation per sub-model at rank one, each stated over every variable of
  # the global model.
  k <- gvar[["model"]][["k"]]
  relations <- length(object[["submodels"]])
  beta <- gvar[["posterior"]][["beta"]][["coeffs"]]

  expect_false(is.null(beta))
  expect_identical(ncol(beta), k * relations)
  expect_identical(gvar[["model"]][["rank"]], relations)
  expect_true(all(is.finite(as.matrix(beta))))

  # The class is prepended, so everything written for the level VAR still
  # dispatches as it did.
  expect_s3_class(gvar, "gvecsolved")
  expect_s3_class(gvar, "bvarmodel")
  expect_identical(attr(beta, "relations"),
                   paste0(names(object[["submodels"]]), "_r1"))
})

test_that("a solved GVAR carries none, and is refused", {
  object <- gvar_estimated(iterations = 30, burnin = 10)
  gvar <- suppressMessages(submodels_to_gvar(object))

  expect_null(gvar[["posterior"]][["beta"]])
  expect_false(inherits(gvar, "gvecsolved"))
  # The bvartools method for a level VAR says why there is nothing to follow.
  expect_error(persistence_profiles(gvar), "no cointegrating relations")
})

test_that("the profiles start at one and are named for their sub-models", {
  object <- gvec_estimated(iterations = 30, burnin = 10, r = 1)
  gvar <- suppressMessages(submodels_to_gvar(object))
  pp <- persistence_profiles(gvar, n_ahead = 6)

  expect_length(pp, length(object[["submodels"]]))
  expect_named(pp, paste0(names(object[["submodels"]]), "_r1"))
  for (relation in pp) {
    expect_identical(colnames(relation), c("median", "lower", "upper"))
    expect_identical(nrow(relation), 7L)
    expect_equal(unname(relation[1, "median"]), 1)
    expect_true(all(relation >= 0))
  }

  draws <- persistence_profiles(gvar, n_ahead = 4, keep_draws = TRUE)
  expect_identical(ncol(draws[[1]]), 5L)
  expect_equal(unname(draws[[1]][, 1]), rep(1, nrow(draws[[1]])))
})

test_that("the global relation is not the sub-model's partial one", {
  # A sub-model can follow only the part of its relation that it governs, and
  # says so; the global model follows the whole relation. The two are therefore
  # allowed to differ, and the point of the method is that they do.
  object <- gvec_estimated(iterations = 30, burnin = 10, r = 1)
  gvar <- suppressMessages(submodels_to_gvar(object))

  global <- persistence_profiles(gvar, n_ahead = 6)
  expect_warning(partial <- persistence_profiles(object[["submodels"]][[1]][[1]],
                                                 n_ahead = 6),
                 "only")
  expect_identical(nrow(partial[[1]]), nrow(global[[1]]))
})

test_that("the vectors grow with the draws when the export is solved in pieces", {
  # The export path stacks the global model a chunk of draws at a time, and the
  # cointegrating vectors are draws like the coefficients: solved in pieces they
  # have to end up with as many rows. Kept as its own test because the fixtures
  # above solve a model held in memory, which is never chunked, so nothing else
  # here would notice a beta left at the size of the first piece.
  object <- gvec_estimated(iterations = 30, burnin = 10, r = 1)

  folder <- file.path(tempdir(), "bgvars-export-persistence")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE)
  write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
  export <- read_export(folder)

  whole <- suppressMessages(submodels_to_gvar(export, chunk = 100))
  pieces <- suppressMessages(submodels_to_gvar(export, chunk = 7))

  for (solved in list(whole, pieces)) {
    beta <- solved[["posterior"]][["beta"]][["coeffs"]]
    expect_identical(nrow(beta), nrow(solved[["posterior"]][["a"]][["coeffs"]]))
    expect_identical(ncol(beta), solved[["model"]][["k"]] * solved[["model"]][["rank"]])
  }
  expect_equal(unclass(pieces[["posterior"]][["beta"]][["coeffs"]]),
               unclass(whole[["posterior"]][["beta"]][["coeffs"]]))
  expect_identical(attr(pieces[["posterior"]][["beta"]][["coeffs"]], "relations"),
                   attr(whole[["posterior"]][["beta"]][["coeffs"]], "relations"))

  # And the profiles can actually be taken from the chunked solve.
  pp <- persistence_profiles(pieces, n_ahead = 4)
  expect_length(pp, solved[["model"]][["rank"]])
  expect_equal(unname(pp[[1]][1, "median"]), 1)

  unlink(folder, recursive = TRUE)
})
