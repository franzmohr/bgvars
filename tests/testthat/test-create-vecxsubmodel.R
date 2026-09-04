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
