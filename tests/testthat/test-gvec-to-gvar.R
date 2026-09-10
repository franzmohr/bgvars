# Rewriting a GVEC model in levels, so that it can be solved.
#
# The sub-models of a GVEC model are in error correction form, which
# submodels_to_gvar() cannot stack. vec_to_var() rewrites one of them in levels
# and gvec_to_gvar() does it for a whole model.

test_that("vec_to_var rewrites a VECX sub-model as a VARX sub-model", {
  object <- gvec_estimated(iterations = 30, burnin = 10)

  model <- vec_to_var(object[["submodels"]][["US"]][[1]])

  expect_s3_class(model, "varxsubmodel")
  expect_s3_class(model, "bvarmodel")
  expect_equal(model[["model"]][["type"]], "VARX")
  expect_equal(model[["model"]][["endogen"]], c("y", "Dp"))
  expect_equal(model[["model"]][["exogen"]], c("y.s", "Dp.s"))
})

test_that("the converted sub-model describes its own coefficients", {
  # submodels_to_gvar() finds the blocks of a sub-model by counting through
  # these entries, so they have to add up to the number of coefficients that is
  # actually there. Getting one of them wrong would not be an error, it would
  # silently read the wrong block.
  object <- gvec_to_gvar(gvec_estimated(iterations = 30, burnin = 10))

  for (submodel in names(object[["submodels"]])) {

    model <- object[["submodels"]][[submodel]][[1]]
    specification <- model[["model"]]

    expected <-
      specification[["k"]] * specification[["k_endogen"]] * specification[["p_endogen"]] +
      specification[["k"]] * specification[["k_exogen"]] * (specification[["p_exogen"]] + 1) +
      specification[["k"]] * specification[["m_global"]] * (specification[["s_global"]] + 1) +
      specification[["k"]] * specification[["n"]]

    expect_equal(ncol(model[["posterior"]][["a"]][["coeffs"]]), expected,
                 info = submodel)
  }
})

test_that("gvec_to_gvar returns a 'gvarmodel' and keeps the model around it", {
  object <- gvec_estimated(iterations = 30, burnin = 10)

  result <- gvec_to_gvar(object)

  expect_s3_class(result, "gvarmodel")
  expect_equal(names(result[["submodels"]]), names(object[["submodels"]]))
  expect_s3_class(result[["submodels"]][["US"]], "modellist")

  # The units, their variables and their weights are not what the
  # transformation is about, and have to come through untouched.
  expect_equal(result[["global"]], object[["global"]])
  expect_equal(result[["weights"]], object[["weights"]])
})

test_that("a converted GVEC model can be solved", {
  object <- gvec_to_gvar(gvec_estimated(iterations = 30, burnin = 10))

  gvec <- submodels_to_gvar(object)

  expect_s3_class(gvec, "bvarmodel")
  expect_equal(gvec[["model"]][["k"]], 6L)
  expect_equal(gvec[["model"]][["endogen"]],
               c("US_y", "US_Dp", "JP_y", "JP_Dp", "CA_y", "CA_Dp"))
  expect_false(gvec[["model"]][["structural"]])
  expect_false(any(is.na(gvec[["posterior"]][["a"]][["coeffs"]])))
})

test_that("a restricted trend becomes a deterministic term of the global model", {
  # gvec_estimated() uses const = "unrestricted" and trend = "restricted", so
  # the level representation carries both.
  object <- gvec_to_gvar(gvec_estimated(iterations = 30, burnin = 10))

  gvec <- submodels_to_gvar(object)

  expect_equal(gvec[["model"]][["n"]], 2L)
  expect_equal(gvec[["model"]][["deterministic"]], c("const", "trend"))
  expect_true(all(c("const", "trend") %in%
                    dimnames(gvec[["data"]][["train"]][["x"]])[[2]]))

  # The trend is the one of the global model, restricted to the estimation
  # sample. It counts from the start of the data rather than from the start of
  # that sample, so that sub-models with different lag orders share it.
  trend <- gvec[["data"]][["train"]][["x"]][, "trend"]
  expect_equal(diff(as.numeric(trend)), rep(1, length(trend) - 1))
  expect_equal(as.numeric(trend),
               as.numeric(stats::window(
                 object[["global"]][["deterministic"]][, "trend"],
                 start = stats::tsp(gvec[["data"]][["train"]][["y"]])[1],
                 end = stats::tsp(gvec[["data"]][["train"]][["y"]])[2])))
  expect_equal(as.numeric(gvec[["data"]][["train"]][["x"]][, "const"]),
               rep(1, length(trend)))
})

test_that("generalised impulse responses of a solved GVEC model are the textbook ones", {
  object <- gvec_to_gvar(gvec_estimated(iterations = 30, burnin = 10))
  gvec <- submodels_to_gvar(object)

  n_ahead <- 5L
  result <- irf(gvec, impulse = "US_y", response = "JP_y",
                n_ahead = n_ahead, ci = 0.68, type = "gir", shock = "sd")

  reference <- reference_girf(gvec, impulse = "US_y", response = "JP_y",
                              n_ahead = n_ahead)
  reference <- t(apply(reference, 1, stats::quantile,
                       probs = c(0.16, 0.5, 0.84)))

  expect_equal(as.matrix(result), reference, ignore_attr = TRUE)
})

test_that("gvec_to_gvar validates its argument", {
  expect_error(gvec_to_gvar(gvar_object()),
               "must be of class 'gvecmodel'")

  expect_error(gvec_to_gvar(gvec_object()),
               "does not contain sub-models")

  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  expect_error(gvec_to_gvar(object), "has not been estimated yet")
})

test_that("the global model spans the sample the sub-models were estimated on", {
  # align_model_obs() trims every sub-model to the longest candidate lag, so a
  # selection of shorter models leaves the global model with room for one more
  # observation than the sub-models actually used. It has to describe their
  # sample, not that longer one.
  object <- gvar_estimated_grid(iterations = 20, burnin = 10)
  object[["submodels"]] <- lapply(object[["submodels"]], function(x) {
    result <- x[1]
    class(result) <- class(x)
    result
  })

  gvar <- submodels_to_gvar(object)

  train <- object[["submodels"]][["US"]][[1]][["data"]][["train"]][["y"]]
  expect_equal(nrow(gvar[["data"]][["train"]][["y"]]), nrow(train))
  expect_equal(stats::tsp(gvar[["data"]][["train"]][["y"]]), stats::tsp(train))
})
