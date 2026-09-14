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
  # p_exogen = 0 and s = 0 put the weakly exogenous and the global variables into
  # the error correction term alone. In levels they enter at lag one with their
  # columns of Pi and at lag zero with a coefficient of zero, which is two
  # blocks. The rewrite counted one and failed on its own column names.
  specs <- list(list(p_exogen = 2, global = "poil", s = 1),
                list(p_exogen = 0, global = NULL, s = NULL),
                list(p_exogen = 1, global = "poil", s = 0),
                list(p_exogen = 0, global = "poil", s = 0))

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

      if (spec$p_exogen == 0) {
        expect_equal(varx[["model"]][["p_exogen"]], 2L, info = info)
        expect_true(all(c("y.s.l00", "y.s.l01") %in% colnames(varx[["data"]][["train"]][["x"]])),
                    info = info)
      }
      if (identical(spec$s, 0)) {
        expect_equal(varx[["model"]][["s_global"]], 2L, info = info)
        expect_true(all(c("poil.l00", "poil.l01") %in% colnames(varx[["data"]][["train"]][["x"]])),
                    info = info)
      }
    }

    gvar <- submodels_to_gvar(object)
    expect_false(anyNA(gvar[["posterior"]][["a"]][["coeffs"]]), info = info)
  }
})
