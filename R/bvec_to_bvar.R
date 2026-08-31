#' Transform a VEC Model to a VAR in Levels
#' 
#' An object of class 'vecxsubmodelest' is transformed to a VAR in levels.
#' 
#' @param object an object of class 'vecxsubmodelest'.
#' 
#' @return An object of class 'varxsubmodelest'.
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
#' @export
bvec_to_bvar.vecxsubmodelest <- function(object) {
  
  # Skip tests if posterior simulation was not successful
  cond <- is.null(object[["error"]])
  
  if (cond) {
    # Get number of draws and draw information
    draws <- NULL
    specs <- NULL
    vars <- c("sigma", "a", "beta")
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
    
    # Model specs  
    k_domestic <- object[["model"]][["k_domestic"]]
    p_domestic <- object[["model"]][["p_domestic"]]
    k_foreign <- object[["model"]][["k_foreign"]]
    p_foreign <- object[["model"]][["p_foreign"]]
    m <- object[["model"]][["m"]]
    s <- object[["model"]][["s"]]
    n_restricted <- object[["model"]][["n_restricted"]]
    n_unrestricted <- object[["model"]][["n_unrestricted"]]
    tt <- nrow(object[["data"]][["y"]])
    tvp <- object[["model"]][["tvp"]]
    r <- object[["model"]][["rank"]]
    n_alpha <- k_domestic * r
    n_gamma_domestic <- max(k_domestic * k_domestic * (p_domestic - 1), 0)
    n_gamma_foreign <- k_domestic * (k_foreign * p_foreign)
    n_upsilon <- k_domestic * (m * s)
    
    # Calculate Pi matrices ----
    object <- .create_pi_matrices(object)
    
    # Produce VAR matrices ----
    
    ## Domestic ----
    domestic <- NULL
    n_a <- k_domestic * k_domestic * p_domestic
    if (tvp) {
      stop("implement tvp")
    } else {
      pos_a_domestic <- n_alpha + 1:n_gamma_domestic
      pos_pi_domestic <- 1:(k_domestic * k_domestic)
    }
    if (p_domestic > 1) {
      
      temp <- bvartools::vec_to_var_transformation_matrix(k_domestic, p_domestic)
      W <- temp[["w"]]
      J <- temp[["j"]]
      
      if (tvp) {
        stop("Implement TVP")
        A <- list()
        for (i in 1:tt) {
          pi_temp <- matrix(0, k, k)
          A[[i]] <- matrix(NA, n_gamma, draws)
          for (draw in 1:draws) {
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
            A[[i]][, draw] <- cbind(pi_temp, gamma_temp) %*% W + J 
          }
          A[[i]] <- t(A[[i]])
        }
        
      } else {
        
        domestic <- matrix(NA, draws, n_a)
        for (draw in 1:draws) {
          if (r > 0) {
            domestic[draw, ] <- cbind(matrix(object[["posteriors"]][["pi"]][["coeffs"]][draw, pos_pi_domestic], k_domestic),
                                      matrix(object[["posteriors"]][["a"]][["coeffs"]][draw, pos_a_domestic], k_domestic)) %*% W + J
          } else {
            domestic[draw, ] <- cbind(matrix(0, k_domestic, k_domestic),
                                      matrix(object[["posteriors"]][["a"]][["coeffs"]][draw, pos_a_domestic], k_domestic)) %*% W + J          
          }
        }
      }
      
    } else {
      
      if (r > 0) {
        if (tvp) {
          stop("implement tvp")
          A <- list()
          for (i in 1:tt) {
            A[[i]] <- matrix(NA, n_a, draws)
            for (draw in 1:draws) {
              A[[i]][, draw] <- matrix(object[["posteriors"]][["pi_domestic"]][[i]][draw, ], k) + matrix(diag(1, k), k)
            }
            A[[i]] <- t(A[[i]])
          }
        } else {
          domestic <- matrix(NA, draws, n_a)
          for (draw in 1:draws) {
            domestic[draw,] <- matrix(object[["posteriors"]][["pi"]][["coeffs"]][draw, pos_pi_domestic], k_domestic) + matrix(diag(1, k_domestic), k_domestic)
          } 
        }
      } else {
        domestic <- matrix(0, draws, n_a)
      }
    }
    
    ## Foreign ----
    foreign <- NULL
    
    W <- bvartools::vec_to_var_transformation_matrix(m = k_foreign, s = p_foreign)[["w_exo"]]
    
    n_b <- k_domestic * k_foreign * (p_foreign + 1)
    if (tvp) {
      stop("implement tvp")
    } else {
      pos_pi_foreign <- k_domestic * k_domestic + 1:(k_domestic * k_foreign)
      pos_a_foreign <- n_alpha + n_gamma_domestic + 1:n_gamma_foreign
    }
    
    if (tvp) {
      stop("implement tvp")
      B <- list()
      for (i in 1:tt) {
        B[[i]] <- matrix(NA, n_b, draws)
        for (draw in 1:draws){
          pix_temp <- matrix(0, k, k_foreign)
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
          B[[i]][, draw] <- cbind(pix_temp, ups_temp) %*% W 
        }
        B[[i]] <- t(B[[i]])
      }
      
    } else {
      foreign <- matrix(NA, draws, n_b)
      for (draw in 1:draws){
        if (r > 0) {
          foreign[draw, ] <- cbind(matrix(object[["posteriors"]][["pi"]][["coeffs"]][draw, pos_pi_foreign], k_domestic),
                                   matrix(object[["posteriors"]][["a"]][["coeffs"]][draw, pos_a_foreign], k_domestic)) %*% W 
        } else {
          foreign[draw, ] <- cbind(matrix(0, k_domestic, k_foreign),
                                   matrix(object[["posteriors"]][["a"]][["coeffs"]][draw, pos_a_foreign], k_domestic)) %*% W
        }
      }
    }
    
    ## Global ----
    global <- NULL
    if (m > 0) {
      
      stop("impelment global variables")
      W <- bvartools::vec_to_var_transformation_matrix(m = m, s = s)[["w_exo"]]
      
      if (tvp) {
        stop("implement tvp")
      } else {
        pos_pi_global <- k_domestic * (k_domestic + k_foreign) + 1:(k_domestic * m)
        pos_a_global <- n_alpha + n_gamma_domestic + n_gamma_foreign + 1:n_upsilon
      }
      
      if (tvp) {
        stop("implement tvp for global")
        B <- list()
        for (i in 1:tt) {
          B[[i]] <- matrix(NA, n_b, draws)
          for (draw in 1:draws){
            pix_temp <- matrix(0, k, m)
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
            B[[i]][, draw] <- cbind(pix_temp, ups_temp) %*% W 
          }
          B[[i]] <- t(B[[i]])
        }
        
      } else {
        global <- matrix(NA, n_upsilon, draws)
        for (draw in 1:draws){
          if (r > 0) {
            global[draw,] <- cbind(matrix(object[["posteriors"]][["pi"]][["coeffs"]][draw, pos_pi_gobal], k_domestic),
                                   matrix(object[["posteriors"]][["a"]][["coeffs"]][draw, pos_a_global], k_domestic)) %*% W 
          } else {
            global[draw,] <- cbind(matrix(0, k_domestic, m),
                                   matrix(object[["posteriors"]][["a"]][["coeffs"]][draw, pos_a_global], k_domestic)) %*% W 
          }
        } 
      }
    }
    
    ## Deterministic ----
    
    deterministic <- NULL
    if (n_restricted + n_unrestricted > 0) {
      
      if (tvp) {
        stop("implement tvp")
      } else {
        pos_pi_det <- k_domestic * (k_domestic + k_foreign + m) + 1:(k_domestic * n_restricted)
        pos_a_det <- n_alpha + n_gamma_domestic + n_gamma_foreign + n_upsilon + 1:(k_domestic * n_unrestricted)
      }
      
      if (tvp) {
        stop("implement tvp")
        object[["posteriors"]][["deterministic"]] <- list()
        for (i in 1:tt) {
          object[["posteriors"]][["deterministic"]][[i]] <- matrix(NA, draws, (k_det_r + k_det_ur) * k)
          if (k_det_r > 0) {
            object[["posteriors"]][["deterministic"]][[i]][, 1:(k_det_r * k)] <- object[["posteriors"]][["pi_deterministic"]][[i]]
          }
          if (k_det_ur > 0) {
            object[["posteriors"]][["deterministic"]][[i]][, (k_det_r * k) + 1:(k_det_ur * k)] <- object[["posteriors"]][["gamma_deterministic"]][[i]]
          }
        }
        object[["posteriors"]][["beta_deterministic"]] <- NULL
        object[["posteriors"]][["gamma_deterministic"]] <- NULL
        
      } else {
        deterministic <- matrix(NA_real_, draws, k_domestic * (n_restricted + n_unrestricted))
        if (n_unrestricted > 0) {
          deterministic[, 1:(n_unrestricted * k_domestic)] <- object[["posteriors"]][["a"]][["coeffs"]][, pos_a_det]
        }
        if (n_restricted > 0 & r > 0) {
          deterministic[, (n_unrestricted * k_domestic) + 1:(n_restricted * k_domestic)] <- object[["posteriors"]][["pi"]][["coeffs"]][, pos_pi_det]
        }
      } 
    }
    
    ## Structural ----
    
    structural <- NULL
    if (object[["model"]][["structural"]]) {
      stop("implement structural")
      if (tvp) {
        stop("implement tvp")
        object[["posteriors"]][["deterministic"]] <- list()
        for (i in 1:tt) {
          object[["posteriors"]][["deterministic"]][[i]] <- matrix(NA, draws, (k_det_r + k_det_ur) * k)
          if (k_det_r > 0) {
            object[["posteriors"]][["deterministic"]][[i]][, 1:(k_det_r * k)] <- object[["posteriors"]][["pi_deterministic"]][[i]]
          }
          if (k_det_ur > 0) {
            object[["posteriors"]][["deterministic"]][[i]][, (k_det_r * k) + 1:(k_det_ur * k)] <- object[["posteriors"]][["gamma_deterministic"]][[i]]
          }
        }
        object[["posteriors"]][["beta_deterministic"]] <- NULL
        object[["posteriors"]][["gamma_deterministic"]] <- NULL
        
      } else {
        deterministic <- matrix(NA_real_, draws, k * (n_restricted + n_unrestricted))
        if (n_unrestricted > 0) {
          pos_a_det <- n_alpha + n_gamma_domestic + n_gamma_foreign + n_upsilon + 1:(k * n_unrestricted)
          deterministic[, 1:(n_unrestricted * k)] <- object[["posteriors"]][["a"]][["coeffs"]][, pos_a_det]
        }
        if (n_restricted > 0 & r > 0) {
          pos_pi_det <- k * (k + k_foreign + m) + 1:(k * n_restricted)
          deterministic[, (n_unrestricted * k) + 1:(n_restricted * k)] <- object[["posteriors"]][["pi"]][["coeffs"]][, pos_pi_det]
        }
      } 
    }
    
    # update a
    object[["posteriors"]][["a"]] <- NULL
    object[["posteriors"]][["beta"]] <- NULL
    object[["posteriors"]][["pi"]] <- NULL
    
    
    if (tvp) {
      stop("implement tvp")
      for (j in 1:tt) {
        object[["posteriors"]][["a"]][["coeffs"]][[j]] <- coda::mcmc(object[["posteriors"]][["a"]][["coeffs"]][[j]])
        attr(object[["posteriors"]][["a"]][["coeffs"]][[j]], "mcpar") <- specs
      }
    } else {
      object[["posteriors"]][["a"]][["coeffs"]] <- cbind(domestic, foreign, global, deterministic, structural)
      object[["posteriors"]][["a"]][["coeffs"]] <- coda::mcmc(object[["posteriors"]][["a"]][["coeffs"]])
      attr(object[["posteriors"]][["a"]][["coeffs"]], "mcpar") <- specs 
    }
    
    
    # Update raw data ----
    
    # Prepare deterministic terms for .gen_varx
    if (n_unrestricted + n_restricted > 0) {
      object[["data"]][["deterministic"]] <- cbind(object[["data"]][["determ_unrestricted"]],
                                                   object[["data"]][["determ_restricted"]])
      det_names <- NULL
      if (n_unrestricted > 0) {
        det_names <- append(det_names, dimnames(object[["data"]][["determ_unrestricted"]])[[2]])  
      }
      if (n_restricted) {
        det_names <- append(det_names, dimnames(object[["data"]][["determ_restricted"]])[[2]])  
      }
      dimnames(object[["data"]][["deterministic"]]) <- list(NULL, det_names)
      
      object[["model"]][["n"]] <- n_unrestricted + n_restricted
      object[["model"]][["deterministic"]] <- det_names
      object[["model"]]["determ_unrestricted"] <- NULL
      object[["model"]]["n_unrestricted"] <- NULL
      object[["model"]]["determ_restricted"] <- NULL
      object[["model"]]["n_restricted"] <- NULL
    }
    
    temp <- .gen_varxsubmodel(object)
    
    object[["data"]][["y"]] <- temp[["y"]]
    object[["data"]][["x"]] <- temp[["x"]]
    object[["data"]][["z"]] <- temp[["z"]]
    object[["data"]][["w"]] <- NULL
    
    object[["model"]][["type"]] <- "VARX"
    
    class(object) <- c("varxsubmodelest", "list") 
  }
  
  return(object) 
}