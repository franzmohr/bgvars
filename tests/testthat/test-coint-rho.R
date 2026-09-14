# The autocorrelation of the state equation of the cointegration space.
#
# Koop, Leon-Gonzalez and Strachan (2011) treat rho as a parameter rather than a
# fixed hyperparameter. The samplers of bvartools draw it as soon as the prior
# on beta carries the support of a uniform prior on it, so what a sub-model
# needs is for add_priors() to check that support and pass it on. Without it
# nothing changes and the posterior carries no rho.

vecx_rho <- function(coint, error = "wishart", iterations = 20, burnin = 10) {
  object <- gvec_object()
  object <- add_submodels(object,
                          endogen = c("y", "Dp"), p_endogen = 2,
                          exogen = c("y", "Dp"), p_exogen = 2,
                          const = "unrestricted", r = 1,
                          tvp = TRUE, error = error,
                          iterations = iterations, burnin = burnin)
  object <- align_model_obs(object)
  sigma <- if (error == "wishart") {
    list(df = 3, scale = 0.0001)
  } else {
    list(shape = 3, rate = 0.01, mu = 0, v_i = 0.01,
         state_variance = 0.05, offset = 1e-8)
  }
  add_priors(object,
             coef = list(v_i = 1, v_i_det = 1 / 10, shape = 3,
                         rate = 1e-8, rate_det = 1e-8),
             coint = coint, sigma = sigma)
}

test_that("the support of rho reaches the prior of every sub-model", {
  drawn <- vecx_rho(list(rho = 0.99, rho_min = 0.9, rho_max = 0.999))
  fixed <- vecx_rho(list(rho = 0.999))

  for (submodel in names(drawn[["submodels"]])) {
    prior <- drawn[["submodels"]][[submodel]][[1]][["priors"]][["beta"]]
    expect_equal(prior[["rho"]], 0.99)
    expect_equal(prior[["rho_min"]], 0.9)
    expect_equal(prior[["rho_max"]], 0.999)

    # Without a support rho stays a fixed hyperparameter, and the sampler must
    # not find half of one.
    prior <- fixed[["submodels"]][[submodel]][[1]][["priors"]][["beta"]]
    expect_null(prior[["rho_min"]])
    expect_null(prior[["rho_max"]])
  }
})

test_that("an incomplete or impossible support of rho is refused", {
  expect_error(vecx_rho(list(rho = 0.99, rho_min = 0.9)),
               "must be specified together")
  expect_error(vecx_rho(list(rho = 0.99, rho_max = 0.999)),
               "must be specified together")
  expect_error(vecx_rho(list(rho = 0.95, rho_min = 0.99, rho_max = 0.9)),
               "0 < rho_min < rho_max <= 1")
  expect_error(vecx_rho(list(rho = 0.95, rho_min = 0.9, rho_max = 1.01)),
               "0 < rho_min < rho_max <= 1")
  expect_error(vecx_rho(list(rho = 0.95, rho_min = 0.99, rho_max = 1)),
               "must lie between")
})

test_that("rho is drawn within its support through the scaled sampler", {
  for (error in c("wishart", "sv+covar")) {

    set.seed(20260912)
    drawn <- vecx_rho(list(rho = 0.99, rho_min = 0.9, rho_max = 0.999), error = error)
    drawn <- add_posterior_coefficients(add_initial_values(drawn))

    set.seed(20260912)
    fixed <- vecx_rho(list(rho = 0.999), error = error)
    fixed <- add_posterior_coefficients(add_initial_values(fixed))

    for (submodel in names(drawn[["submodels"]])) {
      posterior <- drawn[["submodels"]][[submodel]][[1]][["posterior"]][["beta"]]
      rho <- as.numeric(posterior[["rho"]])

      # One draw per retained iteration, inside the support, and moving within
      # it -- a draw that never moves is what a block that silently did nothing
      # would leave behind.
      expect_length(rho, nrow(as.matrix(posterior[["coeffs"]])))
      expect_true(all(rho >= 0.9))
      expect_true(all(rho <= 0.999))
      expect_gt(length(unique(rho)), 1)

      # Held fixed, rho is a hyperparameter the caller already has in the
      # prior, so the posterior hands nothing back.
      expect_null(fixed[["submodels"]][[submodel]][[1]][["posterior"]][["beta"]][["rho"]])
    }
  }
})

test_that("a drawn rho and its support survive an export", {
  set.seed(20260913)
  object <- vecx_rho(list(rho = 0.99, rho_min = 0.9, rho_max = 0.999))
  object <- add_posterior_coefficients(add_initial_values(object))

  folder <- file.path(tempdir(), "bgvars-coint-rho")
  unlink(folder, recursive = TRUE)
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  write_to_hdf5(object, folder = folder)

  original <- object[["submodels"]][["US"]][[1]]
  restored <- read_submodel_from_folder(folder, "US")[[1]]

  expect_equal(as.numeric(restored[["posterior"]][["beta"]][["rho"]]),
               as.numeric(original[["posterior"]][["beta"]][["rho"]]))
  # The reader gives every scalar of a prior back as a 1 x 1 matrix, which the
  # sampler reads as the same number.
  expect_equal(as.numeric(restored[["priors"]][["beta"]][["rho_min"]]), 0.9)
  expect_equal(as.numeric(restored[["priors"]][["beta"]][["rho_max"]]), 0.999)
})
