# The sub-model objects were renamed from the "domestic"/"foreign" wording to
# "endogen"/"exogen" (see the model specifications built by
# create_varxsubmodel and create_vecxsubmodel). Some downstream functions were
# not migrated yet and still read the old names. The tests below describe the
# behaviour that is expected once they are, and are skipped until then so that
# R CMD check stays clean.

test_that("global variables can be used in a GVEC sub-model", {
  skip("create_vecxsubmodel does not build the regressor matrix for global variables")

  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 0,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)

  model <- object[["submodels"]][["US"]][[1]]
  expect_equal(model[["model"]][["m_global"]], 1L)
  expect_true("poil.l00" %in% dimnames(model[["data"]][["train"]][["x"]])[[2]])
  expect_equal(model[["model"]][["m"]],
               model[["model"]][["k_exogen"]] * (model[["model"]][["p_exogen"]] + 1) +
                 model[["model"]][["m_global"]] * (model[["model"]][["s_global"]] + 1))

  expect_no_error(add_priors(object,
                             coef = list(v_i = 0),
                             coint = list(v_i = 0, p_tau_i = 1),
                             sigma = list(df = 3, scale = 0.0001)))
})

test_that("get_submodel_specifications summarises the estimated sub-models", {
  skip("get_submodel_specifications reads the pre-rename model specification")

  object <- gvar_estimated(iterations = 20, burnin = 10, p_endogen = 1:2)

  specs <- get_submodel_specifications(object[["submodels"]][["US"]])

  expect_s3_class(specs, "data.frame")
  expect_equal(nrow(specs), 2L)
  expect_true(all(specs[, "type"] == "VARX"))
  expect_equal(specs[, "lag_domestic"], 1:2)
})

test_that("combine_submodels solves the global model", {
  skip("combine_submodels reads the pre-rename model specification")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  # One model per sub-model is required to solve the global model.
  submodels <- lapply(object[["submodels"]], function(x) x[[1]])
  class(submodels) <- c("submodelestlist", "list")

  gvar <- combine_submodels(submodels)

  expect_s3_class(gvar, "bgvar")
  k <- sum(vapply(submodels, function(x) x[["model"]][["k"]], numeric(1)))
  expect_equal(gvar[["model"]][["k"]], k)
  expect_equal(ncol(gvar[["a0"]]), k^2)
  expect_equal(ncol(gvar[["sigma"]]), k^2)
  expect_equal(nrow(gvar[["a"]]), 20L)
})

test_that("girf computes generalised impulse responses of the global model", {
  skip("requires a working combine_submodels")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  submodels <- lapply(object[["submodels"]], function(x) x[[1]])
  class(submodels) <- c("submodelestlist", "list")
  gvar <- combine_submodels(submodels)

  result <- girf(gvar, impulse = c("US", "y"), response = c("JP", "y"),
                 n.ahead = 5, ci = 0.68)

  expect_s3_class(result, "bgvarirf")
  expect_equal(nrow(result), 6L)
})

test_that("gfevd computes generalised forecast error variance decompositions", {
  skip("requires a working combine_submodels")

  object <- gvar_estimated(iterations = 20, burnin = 10)
  submodels <- lapply(object[["submodels"]], function(x) x[[1]])
  class(submodels) <- c("submodelestlist", "list")
  gvar <- combine_submodels(submodels)

  result <- gfevd(gvar, response = c("US", "y"), n.ahead = 5)

  expect_s3_class(result, "bgvarfevd")
  expect_equal(as.numeric(rowSums(result)), rep(1, nrow(result)))
})

test_that("diff.submodeldata returns an object of class 'submodeldata'", {
  skip("diff.submodeldata drops the class, so the result cannot be passed on")

  submodel_data <- gvar_data()
  result <- diff(submodel_data, variables = "y", multi = 100)

  expect_s3_class(result, "submodeldata")
  expect_no_error(create_gvarmodel(submodel_data = result,
                                   global_data = gvar_global_data()))
})
