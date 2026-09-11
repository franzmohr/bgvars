test_that("create_vecxsubmodel returns one model per specification", {
  object <- gvec_object()

  models <- create_vecxsubmodel(object, submodel = "US",
                                endogen = c("y", "Dp"), p_endogen = 1:2,
                                exogen = c("y", "Dp"), p_exogen = 1,
                                const = "unrestricted", trend = "restricted",
                                r = 0:2,
                                iterations = 10, burnin = 10)

  expect_s3_class(models, "modellist")
  expect_length(models, 2 * 3)
  for (i in models) {
    expect_s3_class(i, "vecxsubmodel")
    expect_s3_class(i, "bvecmodel")
  }

  specs <- t(vapply(models,
                    function(x) unlist(x[["model"]][c("p_endogen", "rank")]),
                    numeric(2)))
  expect_false(any(duplicated(specs)))
  expect_setequal(specs[, "rank"], 0:2)
})

test_that("the model specification reflects the requested variables", {
  object <- gvec_object()
  model <- create_vecxsubmodel(object, submodel = "US",
                               endogen = c("y", "Dp"), p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 1,
                               const = "unrestricted", trend = "restricted",
                               r = 1,
                               iterations = 10, burnin = 10)[[1]]

  specs <- model[["model"]]
  expect_equal(specs[["type"]], "VECX")
  expect_equal(specs[["algorithm"]], "VecNormalWishart")
  expect_equal(specs[["k"]], 2L)
  expect_equal(specs[["k_endogen"]], 2L)
  expect_equal(specs[["p_endogen"]], 1L)
  expect_equal(specs[["k_exogen"]], 2L)
  expect_equal(specs[["p_exogen"]], 1L)
  expect_equal(specs[["rank"]], 1L)
  expect_equal(specs[["endogen"]], c("y", "Dp"))
  expect_equal(specs[["exogen"]], c("y.s", "Dp.s"))
  # The unrestricted constant enters the non-cointegration part, the restricted
  # trend the error correction term.
  expect_equal(specs[["n"]], 1L)
  expect_equal(specs[["n_restricted"]], 1L)
})

test_that("the cointegration space grows with the restricted terms", {
  object <- gvec_object()
  build <- function(...) {
    create_vecxsubmodel(object, submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 1,
                        exogen = c("y", "Dp"), p_exogen = 1, r = 1,
                        iterations = 10, burnin = 10, ...)[[1]]
  }

  # Endogenous variables plus the contemporaneous weakly exogenous variables.
  plain <- build()
  expect_equal(plain[["model"]][["k_beta"]], 4L)
  expect_equal(plain[["model"]][["n_restricted"]], 0L)

  restricted <- build(const = "restricted")
  expect_equal(restricted[["model"]][["k_beta"]], 5L)
  expect_equal(restricted[["model"]][["n_restricted"]], 1L)
  expect_equal(restricted[["model"]][["n"]], 0L)

  unrestricted <- build(const = "unrestricted")
  expect_equal(unrestricted[["model"]][["k_beta"]], 4L)
  expect_equal(unrestricted[["model"]][["n"]], 1L)
})

test_that("data matrices are consistent with the model specification", {
  object <- gvec_object()
  model <- create_vecxsubmodel(object, submodel = "US",
                               endogen = c("y", "Dp"), p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 1,
                               const = "unrestricted", trend = "restricted",
                               r = 1,
                               iterations = 10, burnin = 10)[[1]]

  specs <- model[["model"]]
  y <- model[["data"]][["train"]][["y"]]
  w <- model[["data"]][["train"]][["w"]]
  x <- model[["data"]][["train"]][["x"]]
  z <- model[["data"]][["train"]][["z"]]

  expect_equal(ncol(y), specs[["k"]])
  # The error correction term spans the cointegration space.
  expect_equal(ncol(w), specs[["k_beta"]])
  expect_equal(nrow(w), nrow(y))
  # The SUR representation stacks the non-cointegration regressors and the
  # loadings of the error correction term.
  expect_equal(nrow(z), nrow(y) * specs[["k"]])
  expect_equal(ncol(z), ncol(x) * specs[["k"]] + specs[["k"]] * specs[["rank"]])
})

test_that("create_vecxsubmodel validates its arguments", {
  object <- gvec_object()

  expect_error(create_vecxsubmodel(object, submodel = c("US", "JP")),
               "may only contain one element")
  expect_error(create_vecxsubmodel(object, submodel = "US", global = "poil"),
               "argument 's' must be specified")
  # An unspecified rank makes the function report the ranks it falls back to.
  expect_error(suppressMessages(
    create_vecxsubmodel(object, submodel = "US", error = "unknown")),
    "Invalid specification of argument 'error'")
  expect_error(suppressMessages(
    create_vecxsubmodel(object, submodel = "US", varsel = "unknown")),
    "argument 'varsel' is not supported")
})

test_that("create_vecxsubmodel falls back to all admissible ranks", {
  object <- gvec_object()

  expect_message(
    models <- create_vecxsubmodel(object, submodel = "US",
                                  endogen = c("y", "Dp"), p_endogen = 1,
                                  exogen = c("y", "Dp"), p_exogen = 1,
                                  const = "unrestricted", trend = "restricted",
                                  iterations = 10, burnin = 10),
    "rank")

  ranks <- vapply(models, function(x) x[["model"]][["rank"]], numeric(1))
  expect_equal(sort(unique(ranks)), 0:max(ranks))
})

test_that("add_submodels builds a model list for every sub-model", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted",
                          r = 0:1,
                          iterations = 10, burnin = 10)

  expect_equal(names(object[["submodels"]]),
               unique(object[["global"]][["index"]][, "submodel"]))
  for (i in object[["submodels"]]) {
    expect_s3_class(i, "modellist")
    expect_length(i, 2)
  }
})

test_that("the levels of the global variables enter the error correction term", {
  object <- gvec_object()
  build <- function(...) {
    create_vecxsubmodel(object, submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 1,
                        exogen = c("y", "Dp"), p_exogen = 1, r = 1,
                        const = "unrestricted",
                        iterations = 10, burnin = 10, ...)[[1]]
  }

  # A global block must not change how the other levels are treated: it adds a
  # column to the error correction term rather than replacing it. Counting the
  # levels of the global variables over the ones already there used to be an
  # assignment to n_ect rather than an addition, which left a single column in
  # 'w' and moved every other level into the stationary regressors of 'x'.
  plain <- build()
  with_global <- build(global = "poil", s = 1)

  expect_equal(plain[["model"]][["k_beta"]], 4L)
  expect_equal(dimnames(plain[["data"]][["train"]][["w"]])[[2]],
               c("l.y", "l.Dp", "l.y.s", "l.Dp.s"))

  expect_equal(with_global[["model"]][["k_beta"]], 5L)
  expect_equal(dimnames(with_global[["data"]][["train"]][["w"]])[[2]],
               c("l.y", "l.Dp", "l.y.s", "l.Dp.s", "l.poil"))

  # The non-cointegration term holds differences and the unrestricted constant
  # only. Any level left in there would not be a valid VECX specification.
  x_names <- dimnames(with_global[["data"]][["train"]][["x"]])[[2]]
  expect_equal(x_names, c("d.y.s.00", "d.Dp.s.00", "d.poil.00", "const"))
  expect_equal(ncol(with_global[["data"]][["train"]][["w"]]),
               with_global[["model"]][["k_beta"]])
  expect_equal(ncol(with_global[["data"]][["train"]][["z"]]),
               length(x_names) * with_global[["model"]][["k"]] +
                 with_global[["model"]][["k"]] * with_global[["model"]][["rank"]])
})

test_that("add_submodels keeps the global levels in the error correction term", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 10, burnin = 10)

  # add_submodels.gvecmodel routes through create_vecxsubmodel, so every
  # sub-model must show the same error correction term. The restricted trend
  # adds a sixth column to the five levels.
  for (i in object[["submodels"]]) {
    expect_equal(i[[1]][["model"]][["k_beta"]], 6L)
    expect_equal(dimnames(i[[1]][["data"]][["train"]][["w"]])[[2]],
                 c("l.y", "l.Dp", "l.y.s", "l.Dp.s", "l.poil", "trend"))
  }
})

test_that("a lag order of zero drops the differences but keeps the levels", {
  object <- gvec_object()
  build <- function(...) {
    create_vecxsubmodel(object, submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 1,
                        exogen = c("y", "Dp"), r = 1,
                        const = "unrestricted",
                        iterations = 10, burnin = 10, ...)[[1]]
  }

  # Selecting zero blocks of a regressor used to be written as 1:(n * 0), which
  # is c(1, 0) rather than an empty selection, so the first column of the block
  # was picked up twice and every offset behind it addressed the wrong column.
  # The error correction term is unaffected by the lag orders: the levels of the
  # global variables stay in it whether or not their differences are used.
  no_global_diff <- build(p_exogen = 1, global = "poil", s = 0)
  expect_equal(no_global_diff[["model"]][["k_beta"]], 5L)
  expect_equal(dimnames(no_global_diff[["data"]][["train"]][["w"]])[[2]],
               c("l.y", "l.Dp", "l.y.s", "l.Dp.s", "l.poil"))
  expect_equal(dimnames(no_global_diff[["data"]][["train"]][["x"]])[[2]],
               c("d.y.s.00", "d.Dp.s.00", "const"))
  expect_equal(no_global_diff[["model"]][["m"]], 2L)

  no_exogen_diff <- build(p_exogen = 0, global = "poil", s = 1)
  expect_equal(dimnames(no_exogen_diff[["data"]][["train"]][["x"]])[[2]],
               c("d.poil.00", "const"))
  expect_equal(no_exogen_diff[["model"]][["m"]], 1L)

  # A model whose only regressor is deterministic still has to come together.
  only_det <- build(p_exogen = 0, global = "poil", s = 0)
  expect_equal(dimnames(only_det[["data"]][["train"]][["x"]])[[2]], "const")
  expect_equal(only_det[["model"]][["m"]], 0L)
  expect_equal(only_det[["model"]][["k_beta"]], 5L)

  # Without any regressor at all 'x' is dropped and 'z' holds the loadings of
  # the error correction term only.
  bare <- create_vecxsubmodel(object, submodel = "US",
                              endogen = c("y", "Dp"), p_endogen = 1,
                              exogen = c("y", "Dp"), p_exogen = 0,
                              global = "poil", s = 0, r = 1,
                              iterations = 10, burnin = 10)[[1]]
  expect_null(bare[["data"]][["train"]][["x"]])
  expect_equal(ncol(bare[["data"]][["train"]][["z"]]),
               bare[["model"]][["k"]] * bare[["model"]][["rank"]])
})

test_that("a grid over s keeps the regressors of every model distinct", {
  object <- gvec_object()
  models <- create_vecxsubmodel(object, submodel = "US",
                                endogen = c("y", "Dp"), p_endogen = 1,
                                exogen = c("y", "Dp"), p_exogen = 1,
                                global = "poil", s = 0:2, r = 1,
                                const = "unrestricted",
                                iterations = 10, burnin = 10)

  expect_length(models, 3)
  x_names <- lapply(models, function(i) dimnames(i[["data"]][["train"]][["x"]])[[2]])
  expect_equal(x_names,
               list(c("d.y.s.00", "d.Dp.s.00", "const"),
                    c("d.y.s.00", "d.Dp.s.00", "d.poil.00", "const"),
                    c("d.y.s.00", "d.Dp.s.00", "d.poil.00", "d.poil.01", "const")))

  for (i in models) {
    # A column addressed twice is the signature of the 1:0 selection.
    expect_false(any(duplicated(dimnames(i[["data"]][["train"]][["x"]])[[2]])))
    expect_equal(i[["model"]][["m"]],
                 i[["model"]][["k_exogen"]] * i[["model"]][["p_exogen"]] +
                   i[["model"]][["m_global"]] * i[["model"]][["s_global"]])
    expect_equal(ncol(i[["data"]][["train"]][["z"]]),
                 ncol(i[["data"]][["train"]][["x"]]) * i[["model"]][["k"]] +
                   i[["model"]][["k"]] * i[["model"]][["rank"]])
  }
})

test_that("a sub-model with a global block estimates end to end", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)

  model <- object[["submodels"]][["US"]][[1]]
  expect_equal(model[["model"]][["m_global"]], 1L)
  # A VECX model regresses on differences, so the global variable enters the
  # non-cointegration term as a difference and the error correction term as a
  # level.
  expect_true("d.poil.00" %in% dimnames(model[["data"]][["train"]][["x"]])[[2]])
  expect_true("l.poil" %in% dimnames(model[["data"]][["train"]][["w"]])[[2]])
  # 'p_exogen' and 's_global' count blocks of regressors, the same way they do
  # in a VARX model.
  expect_equal(model[["model"]][["m"]],
               model[["model"]][["k_exogen"]] * model[["model"]][["p_exogen"]] +
                 model[["model"]][["m_global"]] * model[["model"]][["s_global"]])

  expect_no_error(add_priors(object,
                             coef = list(v_i = 0),
                             coint = list(v_i = 0, p_tau_i = 1),
                             sigma = list(df = 3, scale = 0.0001)))
})
