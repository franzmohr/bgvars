test_that("create_varxsubmodel returns one model per lag combination", {
  object <- gvar_object()

  models <- create_varxsubmodel(object, submodel = "US",
                                endogen = c("y", "Dp"), p_endogen = 1:2,
                                exogen = c("y", "Dp"), p_exogen = 0:1,
                                global = "poil", s = 0:1,
                                iterations = 10, burnin = 10)

  expect_s3_class(models, "modellist")
  expect_length(models, 2 * 2 * 2)
  for (i in models) {
    expect_s3_class(i, "varxsubmodel")
    expect_s3_class(i, "bvarmodel")
  }

  specs <- t(vapply(models,
                    function(x) unlist(x[["model"]][c("p_endogen", "p_exogen", "s_global")]),
                    numeric(3)))
  expect_false(any(duplicated(specs)))
  expect_setequal(specs[, "p_endogen"], 1:2)
  expect_setequal(specs[, "s_global"], 0:1)
})

test_that("the model specification reflects the requested variables", {
  object <- gvar_object()
  model <- create_varxsubmodel(object, submodel = "US",
                               endogen = c("y", "Dp"), p_endogen = 2,
                               exogen = c("y", "Dp", "r"), p_exogen = 1,
                               global = "poil", s = 0,
                               deterministic = "const",
                               iterations = 10, burnin = 10)[[1]]

  specs <- model[["model"]]
  expect_equal(specs[["type"]], "VARX")
  expect_equal(specs[["algorithm"]], "VarNormalWishart")
  expect_equal(specs[["k"]], 2L)
  expect_equal(specs[["k_endogen"]], 2L)
  expect_equal(specs[["p_endogen"]], 2L)
  expect_equal(specs[["k_exogen"]], 3L)
  expect_equal(specs[["p_exogen"]], 1L)
  expect_equal(specs[["m_global"]], 1L)
  expect_equal(specs[["s_global"]], 0L)
  expect_equal(specs[["n"]], 1L)
  expect_equal(specs[["endogen"]], c("y", "Dp"))
  expect_equal(specs[["exogen"]], c("y.s", "Dp.s", "r.s"))
  expect_equal(specs[["global"]], "poil")
  expect_false(specs[["structural"]])
  expect_false(specs[["tvp"]])
  # 'm' counts all weakly exogenous and global regressors, lag zero included.
  expect_equal(specs[["m"]],
               specs[["k_exogen"]] * (specs[["p_exogen"]] + 1) +
                 specs[["m_global"]] * (specs[["s_global"]] + 1))
})

test_that("data matrices are consistent with the model specification", {
  object <- gvar_object()
  model <- create_varxsubmodel(object, submodel = "US",
                               endogen = c("y", "Dp"), p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 1,
                               global = "poil", s = 0,
                               deterministic = "const",
                               iterations = 10, burnin = 10)[[1]]

  specs <- model[["model"]]
  y <- model[["data"]][["train"]][["y"]]
  x <- model[["data"]][["train"]][["x"]]
  z <- model[["data"]][["train"]][["z"]]

  expect_equal(ncol(y), specs[["k"]])
  expect_equal(dimnames(y)[[2]], specs[["endogen"]])
  expect_equal(ncol(x),
               specs[["k"]] * specs[["p_endogen"]] + specs[["m"]] + specs[["n"]])
  expect_equal(nrow(x), nrow(y))
  # The SUR representation is the Kronecker product of x and an identity matrix.
  expect_equal(dim(z), c(nrow(x) * specs[["k"]], ncol(x) * specs[["k"]]))
  expect_equal(z, kronecker(x, diag(1, specs[["k"]])), ignore_attr = TRUE)

  # One observation is lost to the lag of the endogenous variables.
  expect_equal(nrow(y), nrow(object[["global"]][["endogen"]]) - 1L)
})

test_that("the endogenous variables are the own series of the sub-model", {
  object <- gvar_object()
  model <- create_varxsubmodel(object, submodel = "JP",
                               endogen = c("y", "Dp"), p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 0,
                               iterations = 10, burnin = 10)[[1]]

  endogen <- model[["data"]][["original"]][["endogen"]]
  expect_equal(dimnames(endogen)[[2]], c("y", "Dp"))
  expect_equal(as.numeric(endogen[, "y"]),
               as.numeric(object[["global"]][["endogen"]][, "JP_y"]))

  # Weakly exogenous variables are marked with a trailing ".s" and are built
  # from foreign series only.
  exogen <- model[["data"]][["original"]][["exogen"]]
  expect_equal(dimnames(exogen)[[2]], c("y.s", "Dp.s"))
  expect_false(any(as.numeric(exogen[, "y.s"]) ==
                     as.numeric(object[["global"]][["endogen"]][, "JP_y"])))
})

test_that("deterministic terms are added as requested", {
  object <- gvar_object()
  build <- function(...) {
    create_varxsubmodel(object, submodel = "US",
                        endogen = c("y", "Dp"), p_endogen = 1,
                        exogen = c("y", "Dp"), p_exogen = 0,
                        iterations = 10, burnin = 10, ...)[[1]]
  }

  none <- build(deterministic = "none")
  expect_equal(none[["model"]][["n"]], 0L)
  expect_null(none[["data"]][["original"]][["deterministic"]])

  const <- build(deterministic = "const")
  expect_equal(dimnames(const[["data"]][["original"]][["deterministic"]])[[2]],
               "const")
  expect_true(all(const[["data"]][["original"]][["deterministic"]][, "const"] == 1))

  both <- build(deterministic = "both")
  expect_equal(dimnames(both[["data"]][["original"]][["deterministic"]])[[2]],
               c("const", "trend"))

  seasonal <- build(deterministic = "const", seasonal = TRUE)
  # Quarterly data produce three seasonal dummies.
  expect_equal(seasonal[["model"]][["n"]], 4L)
  expect_equal(dimnames(seasonal[["data"]][["original"]][["deterministic"]])[[2]],
               c("const", "season.1", "season.2", "season.3"))
})

test_that("a single endogenous variable yields an ARX model", {
  object <- gvar_object()
  model <- create_varxsubmodel(object, submodel = "US",
                               endogen = "y", p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 0,
                               iterations = 10, burnin = 10)[[1]]

  expect_equal(model[["model"]][["type"]], "ARX")
  expect_equal(model[["model"]][["k"]], 1L)
  expect_equal(ncol(model[["data"]][["train"]][["y"]]), 1L)
})

test_that("structural models are flagged and get an additional data block", {
  object <- gvar_object()
  model <- create_varxsubmodel(object, submodel = "US",
                               endogen = c("y", "Dp"), p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 0,
                               structural = TRUE, error = "gamma",
                               iterations = 10, burnin = 10)[[1]]

  expect_equal(model[["model"]][["type"]], "SVARX")
  expect_true(model[["model"]][["structural"]])
  # The structural block adds k * (k - 1) / 2 columns to the SUR matrix.
  k <- model[["model"]][["k"]]
  expect_equal(ncol(model[["data"]][["train"]][["z"]]),
               ncol(model[["data"]][["train"]][["x"]]) * k + k * (k - 1) / 2)
})

test_that("create_varxsubmodel validates its arguments", {
  object <- gvar_object()

  expect_error(create_varxsubmodel(object, submodel = c("US", "JP")),
               "may only contain one element")
  expect_error(create_varxsubmodel(object, submodel = "US", endogen = "unknown"),
               "no variable from argument 'endogen' is available")
  expect_error(create_varxsubmodel(object, submodel = "US", global = "poil"),
               "argument 's' must be specified")
  expect_error(create_varxsubmodel(object, submodel = "US", error = "unknown"),
               "Invalid specification of argument 'error'")
  expect_error(create_varxsubmodel(object, submodel = "US", varsel = "unknown"),
               "argument 'varsel' is not supported")
  expect_error(create_varxsubmodel(object, submodel = "US",
                                   structural = TRUE, error = "wishart"),
               "Structural models cannot be estimated")
})

test_that("add_submodels builds a model list for every sub-model", {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1:2,
                          exogen = c("y", "Dp"), p_exogen = 0,
                          global = "poil", s = 0,
                          iterations = 10, burnin = 10)

  expect_equal(names(object[["submodels"]]),
               unique(object[["global"]][["index"]][, "submodel"]))
  for (i in object[["submodels"]]) {
    expect_s3_class(i, "modellist")
    expect_length(i, 2)
  }

  # Sub-models are built from the same arguments, so their specifications agree.
  expect_equal(object[["submodels"]][["US"]][[1]][["model"]],
               object[["submodels"]][["JP"]][[1]][["model"]])
})

test_that("endogenous variables are used if they are available", {
  object <- gvar_object()

  # A specification is written for every sub-model at once, so it may well name
  # a variable a particular one does not have. That one is dropped rather than
  # objected to, and only a specification that leaves nothing is an error.
  model <- create_varxsubmodel(object, submodel = "US",
                               endogen = c("y", "not_a_variable", "Dp"),
                               p_endogen = 1)[[1]]

  expect_equal(model[["model"]][["endogen"]], c("y", "Dp"))
  expect_equal(model[["model"]][["k_endogen"]], 2L)
})

test_that("the endogenous variables keep the order they were named in", {
  object <- gvar_object()

  model <- create_varxsubmodel(object, submodel = "US",
                               endogen = c("Dp", "y"), p_endogen = 1)[[1]]

  expect_equal(model[["model"]][["endogen"]], c("Dp", "y"))
  expect_equal(dimnames(model[["data"]][["train"]][["y"]])[[2]], c("Dp", "y"))
})
