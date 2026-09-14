test_that("the Minnesota prior of a VECX sub-model reaches the sampler, block by block", {
  # add_priors() stored NULL as the prior precision, and minnesota_prior() let
  # the prior of the deterministic terms overwrite the blocks of the weakly
  # exogenous and global variables. The variances are recalculated here from
  # the formulas in ?minnesota_prior.vecxsubmodel.
  model <- create_vecxsubmodel(gvec_object(), "JP",
                               endogen = c("y", "Dp"), p_endogen = 2,
                               exogen = c("y", "Dp"), p_exogen = 3,
                               global = "poil", s = 3,
                               const = "unrestricted", trend = "restricted", r = 1,
                               iterations = 10, burnin = 10)[[1]]
  kappa <- list(kappa1 = 2, kappa2 = 0.5, kappa3 = 0.5, kappa4 = 5)
  prior <- add_priors(model, coef = list(minnesota = kappa),
                      coint = list(v_i = 0, p_tau_i = 1),
                      sigma = list(df = 3, scale = 0.0001))

  n_a <- ncol(model[["data"]][["train"]][["z"]])
  expect_equal(dim(prior[["priors"]][["a"]][["v_inv"]]), c(n_a, n_a))

  y <- unclass(model[["data"]][["train"]][["y"]])
  w <- unclass(model[["data"]][["train"]][["w"]])
  x <- unclass(model[["data"]][["train"]][["x"]])
  k <- 2
  r <- 1
  vars <- c("y", "Dp")

  # Residual standard deviation of each variable from a regression on its own
  # lagged level, its own lagged difference and every deterministic term.
  s_endo <- vapply(1:k, function(i) {
    X <- cbind(w[, paste0("l.", vars[i])], x[, paste0("d.", vars[i], ".01")],
               w[, "trend"], x[, "const"])
    e <- qr.resid(qr(X), y[, i])
    sqrt(sum(e^2) / (nrow(X) - ncol(X)))
  }, numeric(1))
  s_x <- function(col) stats::sd(x[, col])

  v <- 1 / diag(prior[["priors"]][["a"]][["v_inv"]])
  expect_equal(v[1:(k * r)], rep(kappa$kappa1, k * r))
  V <- matrix(v[-(1:(k * r))], k)
  colnames(V) <- colnames(x)
  # Plain numbers, so that the comparisons below are about values rather than
  # about the names a subset of V carries.
  cell <- function(l, name) unname(V[l, name])

  for (l in 1:k) {
    for (j in 1:k) {
      expected <- if (l == j) {
        kappa$kappa1
      } else {
        kappa$kappa1 * kappa$kappa2 * s_endo[l]^2 / s_endo[j]^2
      }
      expect_equal(cell(l, paste0("d.", vars[j], ".01")), unname(expected))
    }
    for (v_name in c("d.y.s", "d.Dp.s", "d.poil")) {
      base <- unname(kappa$kappa1 * kappa$kappa3 * s_endo[l]^2 / s_x(paste0(v_name, ".00"))^2)
      expect_equal(cell(l, paste0(v_name, ".00")), base, info = v_name)
      expect_equal(cell(l, paste0(v_name, ".01")), base / 4, info = v_name)
    }
    expect_equal(cell(l, "const"), unname(kappa$kappa1 * kappa$kappa4 * s_endo[l]^2))
  }
})
