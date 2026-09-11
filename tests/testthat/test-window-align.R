test_that("window.gvarmodel restricts every sub-model to the given period", {
  object <- gvar_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 1,
                          iterations = 10, burnin = 10)

  result <- stats::window(object, start = c(1990, 1), end = c(2000, 4))

  expect_s3_class(result, "gvarmodel")
  for (i in names(result[["submodels"]])) {
    for (j in seq_along(result[["submodels"]][[i]])) {
      y <- result[["submodels"]][[i]][[j]][["data"]][["train"]][["y"]]
      expect_equal(stats::tsp(y), c(1990, 2000.75, 4))
      expect_equal(nrow(y), 44L)
      expect_equal(nrow(result[["submodels"]][[i]][[j]][["data"]][["train"]][["x"]]),
                   nrow(y))
    }
  }
})

test_that("window.gvecmodel restricts every sub-model to the given period", {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 10, burnin = 10)

  result <- stats::window(object, start = c(1990, 1), end = c(2000, 4))

  expect_s3_class(result, "gvecmodel")
  for (i in names(result[["submodels"]])) {
    y <- result[["submodels"]][[i]][[1]][["data"]][["train"]][["y"]]
    expect_equal(stats::tsp(y), c(1990, 2000.75, 4))
  }
})

test_that("align_model_obs equalises the sample across sub-models", {
  object <- gvar_object()
  # Sub-model specific lag orders leave the sub-models with different samples.
  for (i in c("US", "JP", "CA")) {
    object[["submodels"]][[i]] <-
      create_varxsubmodel(object, submodel = i,
                          endogen = c("y", "Dp"),
                          p_endogen = switch(i, "US" = 1, "JP" = 2, "CA" = 3),
                          exogen = c("y", "Dp"), p_exogen = 1,
                          global = "poil", s = 1,
                          iterations = 10, burnin = 10)
  }

  samples <- function(x) {
    unique(do.call("rbind", lapply(unlist(x[["submodels"]], recursive = FALSE),
                                   function(i) stats::tsp(i[["data"]][["train"]][["y"]]))))
  }

  # Different lag orders consume a different number of initial observations.
  expect_gt(nrow(samples(object)), 1L)

  result <- align_model_obs(object)

  expect_s3_class(result, "gvarmodel")
  aligned <- samples(result)
  expect_equal(nrow(aligned), 1L)
  # The common sample starts where the most heavily lagged model starts and
  # ends where the shortest model ends.
  expect_equal(unname(aligned[1, 1]), max(samples(object)[, 1]))
  expect_equal(unname(aligned[1, 2]), min(samples(object)[, 2]))
})

test_that("align_model_obs equalises the sample of a GVEC model", {
  object <- gvec_object()
  for (i in c("US", "JP", "CA")) {
    object[["submodels"]][[i]] <-
      create_vecxsubmodel(object, submodel = i,
                          endogen = c("y", "Dp"),
                          p_endogen = switch(i, "US" = 1, "JP" = 2, "CA" = 3),
                          exogen = c("y", "Dp"), p_exogen = 1,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 10, burnin = 10)
  }

  result <- align_model_obs(object)

  expect_s3_class(result, "gvecmodel")
  observations <- unlist(lapply(unlist(result[["submodels"]], recursive = FALSE),
                                function(i) nrow(i[["data"]][["train"]][["y"]])))
  expect_length(unique(observations), 1L)
})
