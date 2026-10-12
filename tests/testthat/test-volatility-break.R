# The pandemic specifications a sub-model can carry: the volatility break of
# Lenza and Primiceri (2022) and, for a VARX sub-model, Student t errors. Both
# are BayesTS features of the two constant coefficient Wishart samplers, and
# bvartools checks them; what is tested here is that a sub-model offers them,
# refuses what those samplers refuse, records them in the manifest, and that
# the date a global model is solved at decides which error covariance the
# solution is built from.
#
# What the scales are worth is checked where they are drawn. Here the draws are
# written by hand, so that the arithmetic of the solve is tested and not the
# sampler.

break_object <- function(form = c("vec", "var"), ...) {
  form <- match.arg(form)
  object <- if (form == "vec") gvec_object() else gvar_object()
  args <- list(object, endogen = c("y", "Dp"), p_endogen = 1,
               exogen = c("y", "Dp"), p_exogen = 2,
               iterations = 10, burnin = 5, ...)
  if (form == "vec") {
    args <- c(args, list(r = 1, const = "unrestricted", trend = "restricted"))
  } else {
    args <- c(args, list(deterministic = "const"))
  }
  do.call(add_submodels, args)
}

# 1979.5 is the first period of the estimation sample of the fixtures, which
# run to 2003.75, so this is an ordinary date inside it.
BREAK <- c(1990, 1)

test_that("a sub-model carries a volatility break as a period of its sample", {

  object <- break_object("vec", error = "wishart", volatility_break = BREAK)
  model <- object[["submodels"]][["US"]][[1]]
  y <- model[["data"]][["train"]][["y"]]

  expect_equal(model[["model"]][["volatility_break"]],
               which(abs(stats::time(y) - 1990) < 1e-8))
  expect_equal(model[["model"]][["volatility_scales"]], 3L)

  plain <- break_object("vec", error = "wishart")
  expect_null(plain[["submodels"]][["US"]][[1]][["model"]][["volatility_break"]])
})

test_that("a VARX sub-model carries Student t errors and a VECX refuses them", {

  object <- break_object("var", error = "wishart", student_t = TRUE)
  expect_true(isTRUE(object[["submodels"]][["US"]][[1]][["model"]][["student_t"]]))

  # Not a gap: BayesTS draws the scale of a period in VarNormalWishart and in
  # no error correction sampler, and the argument used to fall into '...' and
  # be dropped without a word.
  expect_error(break_object("vec", error = "wishart", student_t = TRUE),
               "not available for an error correction sub-model")
})

test_that("the break is refused where the sampler could not read it", {

  expect_error(break_object("vec", error = "sv+covar", volatility_break = BREAK),
               "constant coefficients")
  expect_error(break_object("vec", error = "wishart", tvp = TRUE,
                            volatility_break = BREAK),
               "constant coefficients")
  expect_error(break_object("vec", error = "wishart",
                            volatility_break = c(1950, 1)),
               "not a period of the estimation sample")
  expect_error(break_object("var", error = "wishart", student_t = TRUE,
                            tvp = TRUE),
               "constant coefficients")
})

test_that("the prior of a break is required, and refused without one", {

  object <- break_object("vec", error = "wishart", volatility_break = BREAK)
  coef <- list(v_i = 0)
  coint <- list(v_i = 0, p_tau_i = 1)
  expect_error(add_priors(object, coef = coef, coint = coint,
                          sigma = list(df = 3, scale = 1e-4)),
               "volatility")

  with_prior <- add_priors(object, coef = coef, coint = coint,
                           sigma = list(df = 3, scale = 1e-4,
                                        volatility = list(shape = 1, rho_shape1 = 3.04,
                                                          rho_shape2 = 1.51)))
  prior <- with_prior[["submodels"]][["US"]][[1]][["priors"]][["volatility"]]
  expect_equal(sort(names(prior)), c("rho_shape1", "rho_shape2", "shape"))

  plain <- break_object("vec", error = "wishart")
  expect_error(add_priors(plain, coef = coef, coint = coint,
                          sigma = list(df = 3, scale = 1e-4,
                                       volatility = list(shape = 1, rho_shape1 = 3.04,
                                                         rho_shape2 = 1.51))),
               "this model has none")
})

test_that("the date a global model is solved at decides the error covariance", {

  # A sub-model with one covariance and a break: the scales are written by
  # hand, two draws, so that the factor the solve applies is known exactly.
  object <- break_object("vec", error = "wishart", volatility_break = BREAK)
  model <- object[["submodels"]][["US"]][[1]]
  k <- model[["model"]][["k"]]
  first <- model[["model"]][["volatility_break"]]
  draws <- 2L

  model[["posterior"]][["u_sigma_inv"]] <- list(
    "coeffs" = matrix(rep(as.numeric(diag(k)), each = draws), draws))
  # s1 = 2, s2 = s3 = 1, rho = 0.5
  model[["posterior"]][["volatility"]] <- list(
    "coeffs" = matrix(rep(c(2, 1, 1, 0.5), each = draws), draws))

  at <- function(pos) {
    out <- bgvars:::.submodel_at_period(model, stats::time(model[["data"]][["train"]][["y"]])[pos],
                               "US")
    matrix(out[["posterior"]][["u_sigma_inv"]][["coeffs"]][1, ], k)
  }

  # Before the break the covariance is the one that was drawn.
  expect_equal(at(first - 1), diag(k))
  # At the break the covariance is s^2 times it, so the precision is divided
  # by four.
  expect_equal(at(first), diag(k) / 4)
  # The second and third periods have a scale of one again.
  expect_equal(at(first + 1), diag(k))
  # After the free scales the excess decays: s = 1 + (1 - 1) * rho^j = 1 here,
  # the last free scale being one.
  expect_equal(at(first + 3), diag(k))

  # And the solved sub-model no longer carries the break, so nothing can scale
  # it a second time.
  out <- bgvars:::.submodel_at_period(model,
                             stats::time(model[["data"]][["train"]][["y"]])[first], "US")
  expect_null(out[["model"]][["volatility_break"]])
  expect_null(out[["posterior"]][["volatility"]])
})

test_that("the export records what a sub-model carries", {

  object <- break_object("vec", error = "wishart", volatility_break = BREAK)
  object <- add_priors(object, coef = list(v_i = 0), coint = list(v_i = 0, p_tau_i = 1),
                       sigma = list(df = 3, scale = 1e-4,
                                    volatility = list(shape = 1, rho_shape1 = 3.04,
                                                      rho_shape2 = 1.51)))
  object <- add_initial_values(object)
  folder <- file.path(tempdir(), "break-export")
  unlink(folder, recursive = TRUE)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)

  write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
  manifest <- read_export(folder)[["manifest"]]
  expect_true(all(c("volatility_break", "volatility_scales", "student_t") %in%
                    colnames(manifest)))
  expect_equal(as.integer(manifest[1, "volatility_break"]),
               object[["submodels"]][["US"]][[1]][["model"]][["volatility_break"]])
  expect_equal(as.integer(manifest[1, "volatility_scales"]), 3L)
})
