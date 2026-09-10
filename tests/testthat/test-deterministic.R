# Deterministic terms of a global model.
#
# They are built once, by create_gvarmodel()/create_gvecmodel(), on the time
# axis of the global data, and every sub-model takes the ones it uses from
# there. A sub-model that built its own would count a trend from its own first
# observation, and sub-models with different lag orders keep different numbers
# of early observations -- so their trends would be shifted against each other
# and could not be stacked into one global trend.

test_that("the global model provides the deterministic terms", {
  object <- create_gvarmodel(submodel_data = gvar_data(),
                             global_data = gvar_global_data())

  deterministic <- object[["global"]][["deterministic"]]

  expect_s3_class(deterministic, "ts")
  expect_equal(stats::tsp(deterministic),
               stats::tsp(object[["global"]][["endogen"]]))
  expect_equal(dimnames(deterministic)[[2]],
               c("const", "trend", "season.1", "season.2", "season.3"))

  expect_true(all(deterministic[, "const"] == 1))
  expect_equal(as.numeric(deterministic[, "trend"]),
               seq_len(nrow(deterministic)))

  # The seasonal dummies follow the cycle of the data, one per period but the
  # last, and never more than one of them is on.
  seasons <- deterministic[, c("season.1", "season.2", "season.3")]
  expect_true(all(rowSums(seasons) %in% c(0, 1)))
  expect_equal(as.numeric(seasons[stats::cycle(deterministic) == 1, "season.1"]),
               rep(1, sum(stats::cycle(deterministic) == 1)))
})

test_that("a GVEC model provides them as well", {
  object <- create_gvecmodel(submodel_data = gvec_data(),
                             global_data = gvec_global_data())

  expect_equal(dimnames(object[["global"]][["deterministic"]])[[2]],
               c("const", "trend", "season.1", "season.2", "season.3"))
})

test_that("sub-models of different lag order share one trend", {
  # The point of the whole arrangement. Two VARX sub-models of the same unit,
  # one with a lag more than the other, must agree on the trend wherever their
  # samples overlap.
  object <- gvar_object()

  short <- create_varxsubmodel(object, submodel = "US",
                               endogen = c("y", "Dp"), p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 0,
                               deterministic = "both")[[1]]
  long <- create_varxsubmodel(object, submodel = "US",
                              endogen = c("y", "Dp"), p_endogen = 3,
                              exogen = c("y", "Dp"), p_exogen = 0,
                              deterministic = "both")[[1]]

  trend_short <- short[["data"]][["train"]][["x"]][, "trend"]
  trend_long <- long[["data"]][["train"]][["x"]][, "trend"]

  # The longer model loses two more observations at the start, and its trend
  # therefore begins two later -- not at one again.
  expect_equal(as.numeric(utils::head(trend_long, 1)),
               as.numeric(utils::head(trend_short, 1)) + 2)

  overlap <- stats::window(stats::ts(trend_short,
                                     start = stats::tsp(short[["data"]][["train"]][["x"]])[1],
                                     frequency = stats::tsp(short[["data"]][["train"]][["x"]])[3]),
                           start = stats::tsp(long[["data"]][["train"]][["x"]])[1])
  expect_equal(as.numeric(overlap), as.numeric(trend_long))
})

test_that("aligning the sample leaves the sub-models with the same trend", {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1:2,
                          exogen = c("y", "Dp"), p_exogen = 0,
                          deterministic = "both",
                          iterations = 20, burnin = 10)
  object <- align_model_obs(object)

  trends <- unlist(lapply(object[["submodels"]], function(models) {
    lapply(models, function(model) {
      paste0(model[["data"]][["train"]][["x"]][, "trend"], collapse = ",")
    })
  }))

  expect_equal(length(unique(trends)), 1L)
})

test_that("a VECX sub-model takes its trend from the same place", {
  object <- gvec_object()

  model <- create_vecxsubmodel(object, submodel = "US",
                               endogen = c("y", "Dp"), p_endogen = 1,
                               exogen = c("y", "Dp"), p_exogen = 1,
                               const = "unrestricted", trend = "restricted",
                               r = 1)[[1]]

  # The restricted trend sits in the error correction term.
  trend <- model[["data"]][["train"]][["w"]][, "trend"]
  expected <- stats::window(object[["global"]][["deterministic"]][, "trend"],
                            start = stats::tsp(model[["data"]][["train"]][["y"]])[1],
                            end = stats::tsp(model[["data"]][["train"]][["y"]])[2])

  expect_equal(as.numeric(trend), as.numeric(expected))
  expect_true(all(model[["data"]][["train"]][["x"]][, "const"] == 1))
})

test_that("a term the global model does not provide is reported", {
  object <- gvar_object()
  object[["global"]][["deterministic"]] <-
    object[["global"]][["deterministic"]][, "const", drop = FALSE]

  expect_error(create_varxsubmodel(object, submodel = "US",
                                   endogen = c("y", "Dp"), p_endogen = 1,
                                   exogen = c("y", "Dp"), p_exogen = 0,
                                   deterministic = "both"),
               "does not provide the deterministic term")
})
