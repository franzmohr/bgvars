
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
  
  k_dom <- NCOL(object[["data"]][["Y"]])
  k_for <- length(object[["model"]][["foreign"]][["variables"]])
  tt <- NROW(object[["data"]][["Y"]])
  tvp <- object[["model"]][["tvp"]]
  p <- object[["model"]][["domestic"]][["lags"]]
  r <- object[["model"]][["rank"]]
  
  ## Domestic ----
  if (is.null(object[["posteriors"]][["pi_domestic"]])) {
    if (!is.null(object[["posteriors"]][["alpha"]]) & !is.null(object[["posteriors"]][["beta_domestic"]])){
      if (tvp) {
        object[["posteriors"]][["pi_domestic"]] <- list()
        for (j in 1:tt) {
          object[["posteriors"]][["pi_domestic"]][[j]] <- matrix(NA_real_, draws, k_dom * k_dom)
          for (i in 1:draws) {
            object[["posteriors"]][["pi_domestic"]][[j]][i,] <- matrix(object[["posteriors"]][["alpha"]][[j]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_domestic"]][[j]][i,], k_dom)) 
          }
          object[["posteriors"]][["pi_domestic"]][[j]] <- coda::mcmc(object[["posteriors"]][["pi_domestic"]][[j]])
          attr(object[["posteriors"]][["pi_domestic"]][[j]], "mcpar") <- specs
        } 
      } else {
        object[["posteriors"]][["pi_domestic"]] <- matrix(NA_real_, draws, k_dom * k_dom)
        for (i in 1:draws) {
          object[["posteriors"]][["pi_domestic"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_domestic"]][i,], k_dom))
        }
        object[["posteriors"]][["pi_domestic"]] <- coda::mcmc(object[["posteriors"]][["pi_domestic"]])
        attr(object[["posteriors"]][["pi_domestic"]], "mcpar") <- specs
      }
    }
    if (drop_beta) {
      object[["posteriors"]][["beta_domestic"]] <- NULL 
    }
  }
  
  ## Foreign ----
  
  if (is.null(object[["posteriors"]][["pi_foreign"]])) {
    if (!is.null(object[["posteriors"]][["alpha"]]) & !is.null(object[["posteriors"]][["beta_foreign"]])){
      if (tvp) {
        object[["posteriors"]][["pi_foreign"]] <- list()
        for (j in 1:tt) {
          object[["posteriors"]][["pi_foreign"]][[j]] <- matrix(NA_real_, draws, k_dom * k_for)
          for (i in 1:draws) {
            object[["posteriors"]][["pi_foreign"]][[j]][i,] <- matrix(object[["posteriors"]][["alpha"]][[j]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_foreign"]][[j]][i,], k_for)) 
          }
          object[["posteriors"]][["pi_foreign"]][[j]] <- coda::mcmc(object[["posteriors"]][["pi_foreign"]][[j]])
          attr(object[["posteriors"]][["pi_foreign"]][[j]], "mcpar") <- specs
        } 
      } else {
        object[["posteriors"]][["pi_foreign"]] <- matrix(NA_real_, draws, k_dom * k_for)
        for (i in 1:draws) {
          object[["posteriors"]][["pi_foreign"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_foreign"]][i,], k_for))
        }
        object[["posteriors"]][["pi_foreign"]] <- coda::mcmc(object[["posteriors"]][["pi_foreign"]])
        attr(object[["posteriors"]][["pi_foreign"]], "mcpar") <- specs 
      }
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
        if (tvp) {
          object[["posteriors"]][["pi_global"]] <- list()
          for (j in 1:tt) {
            object[["posteriors"]][["pi_global"]][[j]] <- matrix(NA_real_, draws, k_dom * k_glo)
            for (i in 1:draws) {
              object[["posteriors"]][["pi_global"]][[j]][i,] <- matrix(object[["posteriors"]][["alpha"]][[j]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_global"]][[j]][i,], k_glo)) 
            }
            object[["posteriors"]][["pi_global"]][[j]] <- coda::mcmc(object[["posteriors"]][["pi_global"]][[j]])
            attr(object[["posteriors"]][["pi_global"]][[j]], "mcpar") <- specs
          } 
        } else {
          object[["posteriors"]][["pi_global"]] <- matrix(NA_real_, draws, k_dom * k_glo)
          for (i in 1:draws) {
            object[["posteriors"]][["pi_global"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_global"]][i,], k_glo))
          }
          object[["posteriors"]][["pi_global"]] <- coda::mcmc(object[["posteriors"]][["pi_global"]])
          attr(object[["posteriors"]][["pi_global"]], "mcpar") <- specs 
        }
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
        if (tvp) {
          object[["posteriors"]][["pi_deterministic"]] <- list()
          for (j in 1:tt) {
            object[["posteriors"]][["pi_deterministic"]][[j]] <- matrix(NA_real_, draws, k_dom * k_det_r)
            for (i in 1:draws) {
              object[["posteriors"]][["pi_deterministic"]][[j]][i,] <- matrix(object[["posteriors"]][["alpha"]][[j]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_deterministic"]][[j]][i,], k_det_r)) 
            }
            object[["posteriors"]][["pi_deterministic"]][[j]] <- coda::mcmc(object[["posteriors"]][["pi_deterministic"]][[j]])
            attr(object[["posteriors"]][["pi_deterministic"]][[j]], "mcpar") <- specs
          } 
        } else {
          object[["posteriors"]][["pi_deterministic"]] <- matrix(NA_real_, draws, k_dom * k_det_r)
          for (i in 1:draws) {
            object[["posteriors"]][["pi_deterministic"]][i,] <- matrix(object[["posteriors"]][["alpha"]][i,], k_dom) %*% t(matrix(object[["posteriors"]][["beta_deterministic"]][i,], k_det_r))
          }
          object[["posteriors"]][["pi_deterministic"]] <- coda::mcmc(object[["posteriors"]][["pi_deterministic"]])
          attr(object[["posteriors"]][["pi_deterministic"]], "mcpar") <- specs 
        }
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