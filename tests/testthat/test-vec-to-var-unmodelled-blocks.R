# The level representation of a VECX sub-model, block by block.
#
# An independent check of the rewrite: for the same draw, the level model
# y_t = B x_t + u_t and the error correction model
# Dy_t = alpha beta' w_{t-1} + Gamma x_t + u_t must leave the same u_t.

gvec_estimated_blocks <- function(p_exogen, global = NULL, s = NULL) {
  object <- add_submodels(gvec_object(),
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = p_exogen,
                          global = global, s = s,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  object <- align_model_obs(object)
  object <- add_priors(object,
                       coef = list(v_i = 0),
                       coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 0.0001))
  object <- add_initial_values(object)
  add_posterior_coefficients(object)
}

vecx_residuals <- function(model, draw) {
  k <- model[["model"]][["k"]]
  r <- model[["model"]][["rank"]]
  a <- matrix(model[["posterior"]][["a"]][["coeffs"]][draw, ], k)
  beta <- matrix(model[["posterior"]][["beta"]][["coeffs"]][draw, ], ncol = r)
  unclass(model[["data"]][["train"]][["y"]]) -
    unclass(model[["data"]][["train"]][["w"]]) %*% beta %*% t(a[, 1:r, drop = FALSE]) -
    unclass(model[["data"]][["train"]][["x"]]) %*% t(a[, -(1:r), drop = FALSE])
}

varx_residuals <- function(model, draw) {
  k <- model[["model"]][["k"]]
  b <- matrix(model[["posterior"]][["a"]][["coeffs"]][draw, ], k)
  unclass(model[["data"]][["train"]][["y"]]) -
    unclass(model[["data"]][["train"]][["x"]]) %*% t(b)
}

test_that("the level representation has the residuals of the error correction form", {
  # p_exogen = 1 and s = 1 put the weakly exogenous and the global variables into
  # the error correction term alone. In levels they enter at lag one with their
  # columns of Pi and at lag zero with a coefficient of zero, which is two
  # blocks. The rewrite counted one and failed on its own column names.
  specs <- list(list(p_exogen = 3, global = "poil", s = 2),
                list(p_exogen = 1, global = NULL, s = NULL),
                list(p_exogen = 2, global = "poil", s = 1),
                list(p_exogen = 1, global = "poil", s = 1))

  for (spec in specs) {
    info <- paste0("p_exogen = ", spec$p_exogen,
                   ", s = ", if (is.null(spec$s)) "NULL" else spec$s)

    object <- gvec_estimated_blocks(spec$p_exogen, spec$global, spec$s)

    for (i in names(object[["submodels"]])) {
      vecx <- object[["submodels"]][[i]][[1]]
      varx <- vec_to_var(vecx)

      for (draw in c(1, 20)) {
        expect_equal(varx_residuals(varx, draw), vecx_residuals(vecx, draw),
                     ignore_attr = TRUE, info = paste(info, i, draw))
      }

      if (spec$p_exogen == 1) {
        expect_equal(varx[["model"]][["p_exogen"]], 2L, info = info)
        expect_true(all(c("y.s.l00", "y.s.l01") %in% colnames(varx[["data"]][["train"]][["x"]])),
                    info = info)
      }
      if (identical(spec$s, 1)) {
        expect_equal(varx[["model"]][["s_global"]], 2L, info = info)
        expect_true(all(c("poil.l00", "poil.l01") %in% colnames(varx[["data"]][["train"]][["x"]])),
                    info = info)
      }
    }

    gvar <- submodels_to_gvar(object)
    expect_false(anyNA(gvar[["posterior"]][["a"]][["coeffs"]]), info = info)
  }
})

test_that("a structural sub-model keeps its residuals and its A0 in levels", {
  # Regression test. vec_to_var() stopped for structural sub-models. In
  # A0 Dy_t = alpha beta' w_{t-1} + Gamma x_t + u_t the levels satisfy
  # A0 y_t = B x_t + u_t, with A0 in place of the identity on the first lag, and
  # the free elements of A0 stay the last coefficients.
  object <- add_submodels(gvec_object(),
                          endogen = c("y", "Dp", "rs"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", trend = "restricted", r = 1,
                          iterations = 20, burnin = 10)
  object[["submodels"]][["US"]] <- create_vecxsubmodel(object, submodel = "US",
    endogen = c("y", "Dp", "rs"), p_endogen = 2,
    exogen = c("y", "Dp"), p_exogen = 2,
    const = "unrestricted", trend = "restricted", r = 1,
    structural = TRUE, error = "gamma",
    iterations = 20, burnin = 10)
  object <- align_model_obs(object)
  for (i in names(object[["submodels"]])) {
    sigma <- if (i == "US") list(shape = 3, rate = 0.0001) else list(df = 3, scale = 0.0001)
    object[["submodels"]][[i]] <- add_priors(object[["submodels"]][[i]],
                                             coef = list(v_i = 0),
                                             coint = list(v_i = 0, p_tau_i = 1),
                                             sigma = sigma)
  }
  object <- add_initial_values(object)
  object <- add_posterior_coefficients(object)

  vecx <- object[["submodels"]][["US"]][[1]]
  varx <- vec_to_var(vecx)
  k <- vecx[["model"]][["k"]]
  r <- vecx[["model"]][["rank"]]
  n_struct <- k * (k - 1) / 2
  coeffs <- vecx[["posterior"]][["a"]][["coeffs"]]
  n_a <- ncol(coeffs) - n_struct

  expect_true(varx[["model"]][["structural"]])
  for (draw in c(1, 20)) {
    a0 <- diag(1, k)
    a0[lower.tri(a0)] <- coeffs[draw, n_a + 1:n_struct]
    a <- matrix(coeffs[draw, 1:n_a], k)
    beta <- matrix(vecx[["posterior"]][["beta"]][["coeffs"]][draw, ], ncol = r)
    u_vecx <- unclass(vecx[["data"]][["train"]][["y"]]) %*% t(a0) -
      unclass(vecx[["data"]][["train"]][["w"]]) %*% beta %*% t(a[, 1:r, drop = FALSE]) -
      unclass(vecx[["data"]][["train"]][["x"]]) %*% t(a[, -(1:r), drop = FALSE])

    b_all <- varx[["posterior"]][["a"]][["coeffs"]][draw, ]
    b <- matrix(b_all[seq_len(length(b_all) - n_struct)], k)
    u_varx <- unclass(varx[["data"]][["train"]][["y"]]) %*% t(a0) -
      unclass(varx[["data"]][["train"]][["x"]]) %*% t(b)

    expect_equal(u_varx, u_vecx, ignore_attr = TRUE, info = paste("draw", draw))
    expect_equal(b_all[length(b_all) - n_struct + 1:n_struct],
                 coeffs[draw, n_a + 1:n_struct], ignore_attr = TRUE)
  }

  gvec <- submodels_to_gvar(object)
  pos <- match(paste0("US_", vecx[["model"]][["endogen"]]), gvec[["model"]][["endogen"]])
  a0 <- diag(1, k)
  a0[lower.tri(a0)] <- coeffs[1, n_a + 1:n_struct]
  expect_equal(matrix(gvec[["g"]][1, ], gvec[["model"]][["k"]])[pos, pos], a0,
               ignore_attr = TRUE)
})
