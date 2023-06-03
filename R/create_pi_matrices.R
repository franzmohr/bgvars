
# Helper function to obtain Pi matrices from draws of alpha and beta_*

.create_pi_matrices <- function(object, drop_beta = TRUE) {
  
  draws <- NULL
  specs <- NULL
  vars <- c("a0", "alpha",
            "beta_domestic", "beta_foreign", "beta_global", "beta_deterministic",
            "gamma_domestic", "gamma_foreign", "gamma_global", "gamma_deterministic")
  for (i in vars) {
    if (is.null(draws)) {
      if (!is.null(object[["posteriors"]][[i]])) {
        if (is.list(object[["posteriors"]][[i]])) {
          draws <- nrow(object[["posteriors"]][[i]][[1]])
        } else {
          draws <- nrow(object[["posteriors"]][[i]]) 
        }
      }
    }
    if (is.null(specs)) {
      if (is.list(object[["posteriors"]][[i]])) {
        specs <- attr(object[["posteriors"]][[i]][[1]], "mcpar")
      } else {
        specs <- attr(object[["posteriors"]][[i]], "mcpar")
      }
    }
  }
  
  k <- NCOL(object[["data"]][["Y"]])
  k_for <- length(object[["model"]][["foreign"]][["variables"]])
  tt <- NROW(object[["data"]][["Y"]])
  tvp <- object[["model"]][["tvp"]]
  p <- object[["model"]][["domestic"]][["lags"]]
  r <- object[["model"]][["rank"]]

  # use_dom <- !is.null(object[["posteriors"]][["beta_domestic"]])
  # use_for <- !is.null(object[["posteriors"]][["beta_foreign"]])
  # use_glo <- !is.null(object[["posteriors"]][["beta_global"]])
  # if (use_glo) {
  #   k_glo <- length(object[["model"]][["foreign"]][["variables"]])
  # } else {
  #   k_glo <- 0
  # }
  # use_det <- !is.null(object[["posteriors"]][["beta_deterministic"]])
  # if (use_det) {
  #   k_det <- length(object[["model"]][["deterministic"]][["restricted"]])
  # } else {
  #   k_det <- 0
  # }
  # 
  # pos_dom <- 1:k
  # pos_for <- k + 1:k_for
  # if (use_glo) {
  #   pos_glo <- k + k_for + 1:k_glo
  # }
  # if (use_det) {
  #   pos_det <- k + k_for + k_glo + 1:k_det
  # }
  # 
  # Pi <- matrix(NA, NCOL(object[["data"]][["W"]]) * k, draws)
  # beta <- matrix(NA, NCOL(object[["data"]][["W"]]), r)
  # if (is.null(object[["posteriors"]][["alpha"]])) {
  #   stop("Posterior draws of 'alpha' are missing.")
  # }
  # 
  # for (i in 1:draws) {
  # 
  #   beta[pos_dom,] <- object[["posteriors"]][["beta_domestic"]][i,]
  # 
  #   matrix(object[["posteriors"]][["alpha"]][i,], k) %*% t(matrix(, k))
  # }
  # 
  ## Domestic ----
  if (is.null(object[["posteriors"]][["pi_domestic"]])) {
    if (!is.null(object[["posteriors"]][["alpha"]]) & !is.null(object[["posteriors"]][["beta_domestic"]])){
      object[["posteriors"]][["pi_domestic"]] <- matrix(NA_real_, draws, k * k)
      for (i in 1:draws) {
        object[["posteriors"]][["pi_domestic"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k) %*% t(matrix(object[["posteriors"]][["beta_domestic"]][i,], k))
      }
      object[["posteriors"]][["pi_domestic"]] <- coda::mcmc(object[["posteriors"]][["pi_domestic"]])
      attr(object[["posteriors"]][["pi_domestic"]], "mcpar") <- specs
    }
    if (drop_beta) {
      object[["posteriors"]][["beta_domestic"]] <- NULL 
    }
  }
  
  ## Foreign ----
  
  if (is.null(object[["posteriors"]][["pi_foreign"]])) {
    if (!is.null(object[["posteriors"]][["alpha"]]) & !is.null(object[["posteriors"]][["beta_foreign"]])){
      object[["posteriors"]][["pi_foreign"]] <- matrix(NA_real_, draws, k * k_for)
      for (i in 1:draws) {
        object[["posteriors"]][["pi_foreign"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k) %*% t(matrix(object[["posteriors"]][["beta_foreign"]][i,], k_for))
      }
      object[["posteriors"]][["pi_foreign"]] <- coda::mcmc(object[["posteriors"]][["pi_foreign"]])
      attr(object[["posteriors"]][["pi_foreign"]], "mcpar") <- specs
    }
    if (drop_beta) {
      object[["posteriors"]][["beta_foreign"]] <- NULL 
    }
  }
  
  ## Global ----
  global <- !is.null(object[["model"]][["global"]])
  if (global) {
    k_glo <- length(object[["model"]][["global"]][["variables"]])
    if (is.null(object[["posteriors"]][["pi_global"]])) {
      if (!is.null(object[["posteriors"]][["alpha"]]) & !is.null(object[["posteriors"]][["beta_global"]])){
        object[["posteriors"]][["pi_global"]] <- matrix(NA_real_, draws, k * k_glo)
        for (i in 1:draws) {
          object[["posteriors"]][["pi_global"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k) %*% t(matrix(object[["posteriors"]][["beta_global"]][i,], k_glo))
        }
        object[["posteriors"]][["pi_global"]] <- coda::mcmc(object[["posteriors"]][["pi_global"]])
        attr(object[["posteriors"]][["pi_global"]], "mcpar") <- specs
      }
      if (drop_beta) {
        object[["posteriors"]][["beta_global"]] <- NULL 
      }
    } 
  }
  
  ## Deterministic ----
  k_det_r <- 0
  if (!is.null(object[["model"]][["deterministic"]][["restricted"]])) {
    k_det_r <- length(object[["model"]][["deterministic"]][["restricted"]])
    if (is.null(object[["posteriors"]][["pi_deterministic"]])) {
      if (!is.null(object[["posteriors"]][["alpha"]]) & !is.null(object[["posteriors"]][["beta_deterministic"]])){
        object[["posteriors"]][["pi_deterministic"]] <- matrix(NA_real_, draws, k * k_det_r)
        for (i in 1:draws) {
          object[["posteriors"]][["pi_deterministic"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k) %*% t(matrix(object[["posteriors"]][["beta_deterministic"]][i,], k_det_r))
        }
        object[["posteriors"]][["pi_deterministic"]] <- coda::mcmc(object[["posteriors"]][["pi_deterministic"]])
        attr(object[["posteriors"]][["pi_deterministic"]], "mcpar") <- specs
      }
      if (drop_beta) {
        object[["posteriors"]][["beta_deterministic"]] <- NULL 
      }
    }
  }
  
  if (drop_beta) {
    object[["posteriors"]][["alpha"]] <- NULL
  }
  
  return(object)
}