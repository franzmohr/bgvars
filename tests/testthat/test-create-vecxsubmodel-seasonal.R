test_that("seasonal terms are skipped for annual data, as for VARX sub-models", {
  # The warning said that no dummies were generated, and the function then
  # stopped on the dummies it had not generated.
  submodel_data <- gvec_data()
  for (i in names(submodel_data)) {
    e <- submodel_data[[i]][["endogen"]]
    submodel_data[[i]][["endogen"]] <- stats::ts(e[stats::cycle(e) == 4, , drop = FALSE],
                                                 start = 1979, frequency = 1)
  }
  object <- create_gvecmodel(submodel_data = submodel_data)
  object <- add_weight_matrices(object, submodel_data, period = 1999:2001)

  for (seasonal in c("restricted", "unrestricted")) {
    expect_warning(
      model <- create_vecxsubmodel(object, "US", endogen = c("y", "Dp"),
                                   exogen = c("y", "Dp"), r = 1,
                                   const = "unrestricted", seasonal = seasonal,
                                   iterations = 10, burnin = 10)[[1]],
      "No seasonal dummmies")
    expect_false(any(grepl("season", c(colnames(model[["data"]][["train"]][["w"]]),
                                       colnames(model[["data"]][["train"]][["x"]])))),
                 info = seasonal)
    expect_equal(model[["model"]][["deterministic"]], "const", info = seasonal)
    expect_equal(model[["model"]][["n_restricted"]], 0L, info = seasonal)
  }
})
