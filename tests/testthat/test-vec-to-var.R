# Solving a GVEC model.
#
# The sub-models of a GVEC model are in error correction form. vec_to_var()
# rewrites one of them in levels, and submodels_to_gvar() does so for every
# sub-model in that form before it stacks them.

test_that("vec_to_var rewrites a VECX sub-model as a VARX sub-model", {
  object <- gvec_estimated(iterations = 30, burnin = 10)

  model <- vec_to_var(object[["submodels"]][["US"]][[1]])

  expect_s3_class(model, "varxsubmodel")
  expect_s3_class(model, "bvarmodel")
  expect_equal(model[["model"]][["type"]], "VARX")
  expect_equal(model[["model"]][["endogen"]], c("y", "Dp"))
  expect_equal(model[["model"]][["exogen"]], c("y.s", "Dp.s"))
})

test_that("the converted sub-model describes its own coefficients", {
  # submodels_to_gvar() finds the blocks of a sub-model by counting through
  # these entries, so they have to add up to the number of coefficients that is
  # actually there. Getting one of them wrong would not be an error, it would
  # silently read the wrong block.
  object <- gvec_estimated(iterations = 30, burnin = 10)

  for (submodel in names(object[["submodels"]])) {

    model <- vec_to_var(object[["submodels"]][[submodel]][[1]])
    specification <- model[["model"]]

    # Counted exactly as submodels_to_gvar() counts, with no correction: a
    # block of differences spans one more lag in levels than in differences,
    # and it is vec_to_var() that has to say so in 'p_exogen' and 's_global'.
    expected <-
      specification[["k"]] * specification[["k_endogen"]] * specification[["p_endogen"]] +
      specification[["k"]] * specification[["k_exogen"]] * specification[["p_exogen"]] +
      specification[["k"]] * specification[["m_global"]] * specification[["s_global"]] +
      specification[["k"]] * specification[["n"]]

    expect_equal(ncol(model[["posterior"]][["a"]][["coeffs"]]), expected,
                 info = submodel)
  }
})

test_that("a GVEC model is solved directly", {
  gvec <- submodels_to_gvar(gvec_estimated(iterations = 30, burnin = 10))

  expect_s3_class(gvec, "bvarmodel")
  expect_equal(gvec[["model"]][["k"]], 6L)
  expect_equal(gvec[["model"]][["endogen"]],
               c("US_y", "US_Dp", "JP_y", "JP_Dp", "CA_y", "CA_Dp"))
  expect_false(gvec[["model"]][["structural"]])
  expect_false(any(is.na(gvec[["posterior"]][["a"]][["coeffs"]])))
})

test_that("solving a GVEC model is solving its level representation", {
  object <- gvec_estimated(iterations = 30, burnin = 10)

  levels <- object
  levels[["submodels"]] <- lapply(object[["submodels"]], function(x) {
    x[[1]] <- vec_to_var(x[[1]])
    x
  })
  class(levels) <- c("gvarmodel", "list")

  direct <- submodels_to_gvar(object)
  expect_equal(direct, submodels_to_gvar(levels))

  # A model whose sub-models come in both forms, as a selection among VARX and
  # VECX candidates leaves it, is solved as well, and to the same model.
  mixed <- object
  mixed[["submodels"]][["JP"]] <- levels[["submodels"]][["JP"]]
  expect_equal(submodels_to_gvar(mixed), direct)
})

test_that("a restricted trend becomes a deterministic term of the global model", {
  # gvec_estimated() uses const = "unrestricted" and trend = "restricted", so
  # the level representation carries both.
  object <- gvec_estimated(iterations = 30, burnin = 10)

  gvec <- submodels_to_gvar(object)

  expect_equal(gvec[["model"]][["n"]], 2L)
  expect_equal(gvec[["model"]][["deterministic"]], c("const", "trend"))
  expect_true(all(c("const", "trend") %in%
                    dimnames(gvec[["data"]][["train"]][["x"]])[[2]]))

  # The trend is the one of the global model, restricted to the estimation
  # sample. It counts from the start of the data rather than from the start of
  # that sample, so that sub-models with different lag orders share it.
  trend <- gvec[["data"]][["train"]][["x"]][, "trend"]
  expect_equal(diff(as.numeric(trend)), rep(1, length(trend) - 1))
  expect_equal(as.numeric(trend),
               as.numeric(stats::window(
                 object[["global"]][["deterministic"]][, "trend"],
                 start = stats::tsp(gvec[["data"]][["train"]][["y"]])[1],
                 end = stats::tsp(gvec[["data"]][["train"]][["y"]])[2])))
  expect_equal(as.numeric(gvec[["data"]][["train"]][["x"]][, "const"]),
               rep(1, length(trend)))
})

test_that("generalised impulse responses of a solved GVEC model are the textbook ones", {
  gvec <- submodels_to_gvar(gvec_estimated(iterations = 30, burnin = 10))

  n_ahead <- 5L
  result <- irf(gvec, impulse = "US_y", response = "JP_y",
                n_ahead = n_ahead, ci = 0.68, type = "gir", shock = "sd")

  reference <- reference_girf(gvec, impulse = "US_y", response = "JP_y",
                              n_ahead = n_ahead)
  reference <- t(apply(reference, 1, stats::quantile,
                       probs = c(0.16, 0.5, 0.84)))

  expect_equal(as.matrix(result), reference, ignore_attr = TRUE)
})

test_that("submodels_to_gvar reports what it cannot solve", {
  expect_error(submodels_to_gvar(gvec_object()), "does not contain sub-models")

  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  expect_error(submodels_to_gvar(object), "has not been estimated yet")
})

test_that("the global model spans the sample the sub-models were estimated on", {
  # align_model_obs() trims every sub-model to the longest candidate lag, so a
  # selection of shorter models leaves the global model with room for one more
  # observation than the sub-models actually used. It has to describe their
  # sample, not that longer one.
  object <- gvar_estimated_grid(iterations = 20, burnin = 10)
  object[["submodels"]] <- lapply(object[["submodels"]], function(x) {
    result <- x[1]
    class(result) <- class(x)
    result
  })

  gvar <- submodels_to_gvar(object)

  train <- object[["submodels"]][["US"]][[1]][["data"]][["train"]][["y"]]
  expect_equal(nrow(gvar[["data"]][["train"]][["y"]]), nrow(train))
  expect_equal(stats::tsp(gvar[["data"]][["train"]][["y"]]), stats::tsp(train))
})

test_that("a sub-model with global variables is rewritten in levels", {
  # A block of global variables enters the error correction term in levels and
  # the rest of the sub-model in differences, exactly as the weakly exogenous
  # block does, so it has a level representation of its own. The two blocks may
  # carry different lag orders, which is why the transformation is done here
  # rather than by the 'bvecmodel' method.
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 3,
                          global = "poil", s = 2,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 30, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  model <- vec_to_var(object[["submodels"]][["US"]][[1]])
  specification <- model[["model"]]

  expect_s3_class(model, "varxsubmodel")
  expect_equal(specification[["global"]], "poil")
  # The level representation carries the lag orders the sub-model was given.
  expect_equal(specification[["p_exogen"]], 3L)
  expect_equal(specification[["s_global"]], 2L)
  expect_equal(ncol(model[["data"]][["train"]][["x"]]),
               specification[["k_endogen"]] * specification[["p_endogen"]] +
                 specification[["k_exogen"]] * specification[["p_exogen"]] +
                 specification[["m_global"]] * specification[["s_global"]] +
                 specification[["n"]])

  # And a GVEC model with a global block is solved.
  gvec <- submodels_to_gvar(object)
  expect_equal(gvec[["model"]][["exogen"]], "poil")
})

test_that("the level representation fits what the error correction form fits", {
  # The two forms describe the same model, so for every draw the fitted value
  # of the level equation is the fitted value of the difference equation plus
  # the first lag of the endogenous variables. This checks the reconstructed
  # regressors and the transformed coefficients at once.
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 3,
                          global = "poil", s = 2,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  for (submodel in names(object[["submodels"]])) {

    vec <- object[["submodels"]][[submodel]][[1]]
    var <- vec_to_var(vec)

    k <- vec[["model"]][["k_endogen"]]
    r <- vec[["model"]][["rank"]]
    w <- as.matrix(vec[["data"]][["train"]][["w"]])
    x_vec <- as.matrix(vec[["data"]][["train"]][["x"]])
    x_var <- as.matrix(var[["data"]][["train"]][["x"]])

    for (draw in seq_len(nrow(vec[["posterior"]][["a"]][["coeffs"]]))) {
      a <- matrix(vec[["posterior"]][["a"]][["coeffs"]][draw, ], k)
      beta <- matrix(vec[["posterior"]][["beta"]][["coeffs"]][draw, ], ncol(w), r)
      pi_plus <- a[, seq_len(r), drop = FALSE] %*% t(beta)

      fitted_vec <- w %*% t(pi_plus) +
        x_vec %*% t(a[, r + seq_len(ncol(x_vec)), drop = FALSE])
      fitted_var <- x_var %*%
        t(matrix(var[["posterior"]][["a"]][["coeffs"]][draw, ], k))

      expect_equal(fitted_var, fitted_vec + w[, 1:k], info = submodel,
                   ignore_attr = TRUE)
    }
  }
})

test_that("sub-models with different deterministic terms can be solved together", {
  # A restricted trend drops out of the level representation of a sub-model at
  # rank zero, so it has one deterministic term fewer than the sub-models of
  # rank one. The global model uses both terms and gives that sub-model a trend
  # coefficient of zero.
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 1,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", trend = "restricted", r = 0:1,
                          iterations = 20, burnin = 10)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  # US at rank zero, JP and CA at rank one.
  rank_wanted <- c("US" = 0, "JP" = 1, "CA" = 1)
  for (submodel in names(object[["submodels"]])) {
    models <- object[["submodels"]][[submodel]]
    ranks <- vapply(models, function(x) as.numeric(x[["model"]][["rank"]]), numeric(1))
    result <- models[ranks == rank_wanted[[submodel]]]
    class(result) <- class(models)
    object[["submodels"]][[submodel]] <- result
  }

  # The level representations the global model is stacked from.
  levels <- lapply(object[["submodels"]], function(x) vec_to_var(x[[1]]))

  expect_equal(levels[["US"]][["model"]][["deterministic"]], "const")
  expect_equal(levels[["JP"]][["model"]][["deterministic"]], c("const", "trend"))

  gvec <- submodels_to_gvar(object)

  expect_equal(gvec[["model"]][["n"]], 2L)
  expect_equal(gvec[["model"]][["deterministic"]], c("const", "trend"))
  expect_equal(utils::tail(dimnames(gvec[["data"]][["train"]][["x"]])[[2]], 2),
               c("const", "trend"))

  k <- gvec[["model"]][["k"]]
  n <- gvec[["model"]][["n"]]
  coeffs <- gvec[["posterior"]][["a"]][["coeffs"]]
  pos_c <- ncol(coeffs) - k * n + 1:(k * n)
  endogen <- gvec[["model"]][["endogen"]]

  for (draw in seq_len(nrow(coeffs))) {
    # Back to the structural form, where each block is one sub-model.
    g <- matrix(gvec[["g"]][draw, ], k)
    c_structural <- g %*% matrix(coeffs[draw, pos_c], k)

    for (submodel in names(levels)) {
      model <- levels[[submodel]]
      k_i <- model[["model"]][["k_endogen"]]
      n_i <- model[["model"]][["n"]]
      a_i <- model[["posterior"]][["a"]][["coeffs"]]
      c_i <- matrix(a_i[draw, ncol(a_i) - k_i * n_i + 1:(k_i * n_i)], k_i)
      dimnames(c_i) <- list(NULL, model[["model"]][["deterministic"]])

      rows <- startsWith(endogen, paste0(submodel, "_"))
      expect_equal(c_structural[rows, 1], c_i[, "const"], info = submodel)
      if (submodel == "US") {
        expect_equal(c_structural[rows, 2], rep(0, k_i), info = submodel)
      } else {
        expect_equal(c_structural[rows, 2], c_i[, "trend"], info = submodel)
      }
    }
  }
})
