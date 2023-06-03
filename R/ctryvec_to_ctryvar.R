#' Transform the BVEC Submodels of a BGVAR model to a BVAR Submodel in Levels
#' 
#' An object of class \code{"bgvecest"} is transformed to a VAR in level representation.
#' 
#' @param object an object of class \code{"bgvecest"}.
#' 
#' @return An object of class \code{"bgvarest"}.
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
ctryvec_to_ctryvar <- function(object) {
  
  if (!any(class(object) %in% "bgvecest")) {
    stop("Argument 'object' must be of class 'bgvecest'.")
  }
  
  object <- lapply(object, .transform_ctryvec)
  
  class(object) <- append("bgvarest", class(object)) 
  
  return(object) 
}

.transform_ctryvec <- function(object) {
  
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
  tt <- NROW(object[["data"]][["Y"]])
  tvp <- object[["model"]][["tvp"]]
  p <- object[["model"]][["domestic"]][["lags"]]
  r <- object[["model"]][["rank"]]
  
  # Calculate Pi matrices ----
  
  object <- .create_pi_matrices(object)
  
  # Produce VAR matrices ----
  
  ## Domestic ----
  A <- NULL
  if (!is.null(object[["posteriors"]][["gamma_domestic"]])) {
    
    W <- diag(-1, k * p)
    W[1:k, 1:k] <- diag(1, k)
    W[-(1:k), -(k * (p - 1) + 1:k)] <- W[-(1:k),-(k * (p - 1) + 1:k)] + diag(k * (p - 1))
    J <- matrix(0, k, k * p)
    J[1:k, 1:k] <- diag(1, k)
    
    n_gamma <- k * k * p
    
    if (!tvp) {
      
      A <- matrix(NA, n_gamma, draws)
      for (draw in 1:draws) {
        if (is.null(object[["posteriors"]][["pi_domestic"]])) {
          A[, draw] <- cbind(matrix(0, k, k), matrix(object[["posteriors"]][["gamma_domestic"]][draw, ], k)) %*% W + J
        } else {
          A[, draw] <- cbind(matrix(object[["posteriors"]][["pi_domestic"]][draw, ], k), matrix(object[["posteriors"]][["gamma_domestic"]][draw, ], k)) %*% W + J 
        }
      } 
      
    } else {
      
      A <- matrix(NA, n_gamma * tt, draws)
      pi_temp <- matrix(0, k, k)
      for (draw in 1:draws) {
        for (i in 1:tt) {
          if (!is.null(object[["posteriors"]][["pi_domestic"]])) {
            if (tvp) {
              pi_temp <- matrix(object[["posteriors"]][["pi_domestic"]][[i]][draw, ], k)
            } else {
              pi_temp <- matrix(object[["posteriors"]][["pi_domestic"]][draw, ], k)
            }
          }
          if (tvp) {
            gamma_temp <- matrix(object[["posteriors"]][["gamma_domestic"]][[i]][draw, ], k)
          } else {
            gamma_temp <- matrix(object[["posteriors"]][["gamma_domestic"]][draw, ], k)
          }
          A[(i - 1) * n_gamma + 1:n_gamma, draw] <- cbind(pi_temp, gamma_temp) %*% W + J 
        }
      }
      
    }
    
  } else {
    
    if (!is.null(object[["posteriors"]][["pi_domestic"]])) {
      n_a <- k * k
      if (tvp) {
        A <- matrix(NA, n_a * tt, draws)
        for (draw in 1:draws) {
          for (i in 1:tt) {
            A[(i - 1) * n_a + 1:n_a, draw] <- matrix(object[["posteriors"]][["pi_domestic"]][[i]][draw, ], k) + matrix(diag(1, k), k)
          }
        }
      } else {
        A <- matrix(NA, n_a, draws)
        for (draw in 1:draws) {
          A[, draw] <- matrix(object[["posteriors"]][["pi_domestic"]][draw, ], k) + matrix(diag(1, k), k)
        } 
      }
    } else {
      A <- matrix(0, k * k, draws)
    }
  }
  
  object[["posteriors"]][["domestic"]] <- t(A)
  rm(A)
  object[["posteriors"]][["pi_domestic"]] <- NULL
  object[["posteriors"]][["gamma_domestic"]] <- NULL
  
  ## Foreign ----
  B <- NULL
  if (!is.null(object[["posteriors"]][["gamma_foreign"]])) {
    
    k_for <- length(object[["model"]][["foreign"]][["variables"]])
    p_for <- object[["model"]][["foreign"]][["lags"]]
    
    W <- diag(-1, k_for * (p_for + 1))
    W[1:k_for, 1:k_for] <- 0
    W[1:k_for, k_for + 1:k_for] <- diag(1, k_for)
    W[-(1:k_for), 1:(k_for * p_for)] <- W[-(1:k_for), 1:(k_for * p_for)] + diag(1, k_for * p_for)
    
    n_b <- k * k_for * (p_for + 1)
    
    if (!tvp) {
      B <- matrix(NA, n_b, draws)
      for (draw in 1:draws){
        if (!is.null(object[["posteriors"]][["pi_foreign"]])) {
          B[, draw] <- cbind(matrix(object[["posteriors"]][["pi_foreign"]][draw, ], k), matrix(object[["posteriors"]][["gamma_foreign"]][draw, ], k)) %*% W 
        } else {
          B[, draw] <- cbind(matrix(0, k, k_for), matrix(object[["posteriors"]][["gamma_foreign"]][draw, ], k)) %*% W 
        }
      } 
      
    } else {
      
      B <- matrix(NA, n_b * tt, draws)
      pix_temp <- matrix(0, k, k_for)
      for (draw in 1:draws){
        for (i in 1:tt) {
          if (!is.null(object[["posteriors"]][["pi_foreign"]])) {
            if (tvp) {
              pix_temp <- matrix(object[["posteriors"]][["pi_foreign"]][[i]][draw, ], k)
            } else {
              pix_temp <- matrix(object[["posteriors"]][["pi_foreign"]][draw, ], k)
            } 
          }
          if (tvp) {
            ups_temp <- matrix(object[["posteriors"]][["gamma_foreign"]][[i]][draw, ], k)
          } else {
            ups_temp <- matrix(object[["posteriors"]][["gamma_foreign"]][draw, ], k)
          }
          B[(i - 1) * n_b + 1:n_b, draw] <- cbind(pix_temp, ups_temp) %*% W 
        }
      }
    }
  }
  
  object[["posteriors"]][["foreign"]] <- t(B)
  rm(B)
  object[["posteriors"]][["pi_foreign"]] <- NULL
  object[["posteriors"]][["gamma_foreign"]] <- NULL
  
  ## Global ----
  B <- NULL
  if (!is.null(object[["posteriors"]][["gamma_global"]])) {
    
    k_glo <- length(object[["model"]][["global"]][["variables"]])
    p_glo <- object[["model"]][["global"]][["lags"]]
    
    W <- diag(-1, k_glo * (p_glo + 1))
    W[1:k_glo, 1:k_glo] <- 0
    W[1:k_glo, k_glo + 1:k_glo] <- diag(1, k_glo)
    W[-(1:k_glo), 1:(k_glo * p_glo)] <- W[-(1:k_glo), 1:(k_glo * p_glo)] + diag(1, k_glo * p_glo)
    
    n_b <- k * k_glo * (p_glo + 1)
    
    if (!tvp) {
      B <- matrix(NA, n_b, draws)
      for (draw in 1:draws){
        if (!is.null(object[["posteriors"]][["pi_global"]])) {
          B[, draw] <- cbind(matrix(object[["posteriors"]][["pi_global"]][draw, ], k), matrix(object[["posteriors"]][["gamma_global"]][draw, ], k)) %*% W 
        } else {
          B[, draw] <- cbind(matrix(0, k, k_for), matrix(object[["posteriors"]][["gamma_global"]][draw, ], k)) %*% W 
        }
      } 
      
    } else {
      
      B <- matrix(NA, n_b * tt, draws)
      pix_temp <- matrix(0, k, k_glo)
      for (draw in 1:draws){
        for (i in 1:tt) {
          if (!is.null(object[["posteriors"]][["pi_global"]])) {
            if (tvp) {
              pix_temp <- matrix(object[["posteriors"]][["pi_global"]][[i]][draw, ], k)
            } else {
              pix_temp <- matrix(object[["posteriors"]][["pi_global"]][draw, ], k)
            } 
          }
          if (tvp) {
            ups_temp <- matrix(object[["posteriors"]][["gamma_global"]][[i]][draw, ], k)
          } else {
            ups_temp <- matrix(object[["posteriors"]][["gamma_global"]][draw, ], k)
          }
          B[(i - 1) * n_b + 1:n_b, draw] <- cbind(pix_temp, ups_temp) %*% W 
        }
      }
    }
    object[["posteriors"]][["global"]] <- t(B)
  }
  rm(B)
  object[["posteriors"]][["pi_global"]] <- NULL
  object[["posteriors"]][["gamma_global"]] <- NULL
  
  ## Deterministic ----
  k_det_r <- 0
  if (!is.null(object[["model"]][["deterministic"]][["restricted"]])) {
    k_det_r <- length(object[["model"]][["deterministic"]][["restricted"]])
  }
  k_det_ur <- 0
  if (!is.null(object[["model"]][["deterministic"]][["unrestricted"]])) {
    k_det_ur <- length(object[["model"]][["deterministic"]][["unrestricted"]])
  }
  
  if (k_det_r + k_det_ur > 0) {
    if (tvp) {
      # To be implemented
    } else {
      object[["posteriors"]][["deterministic"]] <- matrix(NA_real_, draws, (k_det_r + k_det_ur) * k)
      if (k_det_r > 0) {
        object[["posteriors"]][["deterministic"]][, 1:(k_det_r * k)] <- object[["posteriors"]][["pi_deterministic"]]
        object[["posteriors"]][["pi_deterministic"]] <- NULL
      }
      if (k_det_ur > 0) {
        object[["posteriors"]][["deterministic"]][, (k_det_r * k) + 1:(k_det_ur * k)] <- object[["posteriors"]][["gamma_deterministic"]]
        object[["posteriors"]][["gamma_deterministic"]] <- NULL
      }
    } 
  }
  
  for (i in c("a0", "domestic", "foreign", "global", "deterministic")) {
    if (!is.null(object[["posteriors"]][[i]])) {
      object[["posteriors"]][[i]] <- coda::mcmc(object[["posteriors"]][[i]])
      attr(object[["posteriors"]][[i]], "mcpar") <- specs
    }
  }
  
  # Update raw data ----
  
  # Prepare deterministic terms for .gen_varx
  if (k_det_r + k_det_ur > 0) {
    object[["data"]][["deterministic"]] <- cbind(object[["data"]][["deterministic"]][["restricted"]],
                                                 object[["data"]][["deterministic"]][["unrestricted"]])
  }
  temp <- .gen_varx(object)
  
  object[["data"]][["Y"]] <- temp[["Y"]]
  object[["data"]][["Z"]] <- temp[["Z"]]
  object[["data"]][["SUR"]] <- temp[["SUR"]]
  object[["data"]][["W"]] <- NULL
  object[["data"]][["X"]] <- NULL
  
  object[["model"]][["type"]] <- "VAR"
  
  class(object) <- c("ctryvarest", "list")
  
  return(object)
}
