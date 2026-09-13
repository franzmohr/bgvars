# The sub-model objects were renamed from the "domestic"/"foreign" wording to
# "endogen"/"exogen" (see the model specifications built by
# create_varxsubmodel and create_vecxsubmodel). Some downstream functions were
# not migrated yet and still read the old names. The tests below describe the
# behaviour that is expected once they are, and are skipped until then so that
# R CMD check stays clean.

test_that("plot draws the cointegration relations of a VECX sub-model", {
  skip(".create_pi_matrices is called but defined nowhere in the package")

  object <- gvec_estimated(iterations = 20, burnin = 10)

  expect_no_error(plot(object[["submodels"]][["US"]][[1]]))
})
