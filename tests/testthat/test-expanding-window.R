# Expanding window estimation of a sub-model.
#
# A sub-model is a 'bvarmodel' or a 'bvecmodel', so the methods of bvartools
# prepare it. bgvars had a method of its own for VARX sub-models, which read
# fields that do not exist and failed on every call.

test_that("sub-models are prepared for expanding window estimation", {
  models <- list(
    varx = create_varxsubmodel(gvar_object(), "US", endogen = c("y", "Dp"),
                               exogen = c("y", "Dp"), p_exogen = 1,
                               iterations = 10, burnin = 10)[[1]],
    vecx = create_vecxsubmodel(gvec_object(), "US", endogen = c("y", "Dp"),
                               exogen = c("y", "Dp"), r = 1, const = "unrestricted",
                               iterations = 10, burnin = 10)[[1]])

  for (name in names(models)) {
    model <- models[[name]]
    y <- model[["data"]][["train"]][["y"]]
    tt <- nrow(y)
    start <- stats::time(y)[tt - 4]

    windows <- use_expanding_window(model, start = start)

    expect_s3_class(windows, "expandingwindow")
    expect_length(windows, 5L + 1L)
    expect_equal(nrow(windows[[1]][["data"]][["train"]][["y"]]), tt - 5L, info = name)
    expect_equal(nrow(windows[[6]][["data"]][["train"]][["y"]]), tt, info = name)
    expect_equal(nrow(windows[[1]][["data"]][["train"]][["x"]]), tt - 5L, info = name)
    expect_equal(nrow(windows[[1]][["data"]][["train"]][["z"]]), 2L * (tt - 5L), info = name)
    expect_s3_class(windows[[1]], class(model)[1])
  }
})
