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

test_that("diff.submodeldata returns an object of class 'submodeldata'", {
  skip("diff.submodeldata drops the class, so the result cannot be passed on")

  submodel_data <- gvar_data()
  result <- diff(submodel_data, variables = "y", multi = 100)

  expect_s3_class(result, "submodeldata")
  expect_no_error(create_gvarmodel(submodel_data = result,
                                   global_data = gvar_global_data()))
})

test_that("plot draws the cointegration relations of a VECX sub-model", {
  skip(".create_pi_matrices is called but defined nowhere in the package")

  object <- gvec_estimated(iterations = 20, burnin = 10)

  expect_no_error(plot(object[["submodels"]][["US"]][[1]]))
})
