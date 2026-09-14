test_that("print shows the sub-models of an estimated global model", {
  # print.gvarmodel stopped with a note to update it as soon as a model held
  # sub-models, and a 'gvecmodel' had no print method at all.
  gvar <- gvar_estimated(iterations = 20, burnin = 10)
  gvec <- gvec_estimated(iterations = 20, burnin = 10)

  for (object in list(gvar, gvec)) {
    output <- capture.output(result <- withVisible(print(object)))
    expect_false(result$visible)
    expect_identical(result$value, object)
    expect_true(any(grepl("Sub-models:", output, fixed = TRUE)))
    expect_true(any(grepl("Weight matrices: added", output, fixed = TRUE)))
  }

  output <- capture.output(print(gvec))
  expect_true(any(grepl("Global Vector Error Correction Model", output, fixed = TRUE)))
  table <- output[(grep("Sub-models:", output, fixed = TRUE) + 2):length(output)]
  expect_match(table[1], "rank")
  expect_length(table, 1 + length(gvec[["submodels"]]))

  output <- capture.output(print(gvar))
  expect_true(any(grepl("Global Vector Autoregressive Model", output, fixed = TRUE)))
  expect_false(any(grepl("rank", output, fixed = TRUE)))
})
