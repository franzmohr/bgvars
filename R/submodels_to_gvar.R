#' Generate Global Model from Sub-Models
#' 
#' ddadsfsf
#' 
#' @param object
#' 
#' @export
submodels_to_gvar <- function(object) {
  
  index <- object[["global"]][["index"]]
  submodels <- unique(index[, "submodel"])
  
  # Number of models per sub-model
  n_models_i <- unlist(lapply(object[["submodels"]], function(x){length(x)}))
  if (any(n_models_i > 1)) {
    pos <- which(n_models_i > 1)
    msg_names <- paste0(names(n_models_i)[pos], collapse = ", ")
    stop("Argument 'object' may only contain one model per sub-model.\n",
         "Affected sub-models: ", msg_names, ".")
  }
  
  # Endogenous variables
  endogen_i <- lapply(object[["submodels"]], function(x){x[[1]][["model"]][["endogen"]]})
  pos_i <- endogen_i
  for (i in submodels) {
    index_i <- index[index[, "submodel"] == i & index[, "variable"] %in% endogen_i[[i]], ]
    pos_i[[i]] <- index_i[match(index_i[,"variable"], endogen_i[[i]]), "id"]
    rm(index_i)
    endogen_i[[i]] <- paste0(i, "_", endogen_i[[i]])
  }
  
  # Number of endogenous variables
  k_endogen_i <- unlist(lapply(object[["submodels"]],
                               function(x){x[[1]][["model"]][["k_endogen"]]}))
  # Maximum lag of endogenous variables
  p_endogen_i <- unlist(lapply(object[["submodels"]],
                               function(x){x[[1]][["model"]][["p_endogen"]]}))
  n_endogen_i <- k_endogen_i * k_endogen_i * p_endogen_i
  
  # Number of weakly exogenous variables
  k_exogen_i <- unlist(lapply(object[["submodels"]],
                              function(x){x[[1]][["model"]][["k_exogen"]]}))
  # Maximum lag of exogenous variables
  p_exogen_i <- unlist(lapply(object[["submodels"]],
                              function(x){x[[1]][["model"]][["p_exogen"]]}))
  n_exogen_i <- k_endogen_i * k_exogen_i * (p_exogen_i + 1)
  
  # Final dimensions of the global model
  k <- sum(k_endogen_i)
  p <- max(p_endogen_i, p_exogen_i)
  k_pos <- cumsum(k_endogen_i)
  k_pos[2:length(k_pos)] <- k_pos[-length(k_pos)]
  k_pos[1] <- 0
  
  # Use global variables
  m_i <- unlist(lapply(object[["submodels"]],
                       function(x){x[[1]][["model"]][["m_global"]]}))
  m <- unique(m_i)
  if (length(m) > 1) {
    stop("Number of global variables differs across sub-models.\n",
         "Feel free to send a feature request.")
  }
  global <- m > 0
  if (global) {
    s_i <- unlist(lapply(object[["submodels"]],
                         function(x){x[[1]][["model"]][["s_global"]]}))
    s <- max(s_i) # Lag of global variables in the global model
    n_global_i <- k_endogen_i * m_i * (s_i + 1)
    global_i <- unique(unlist(lapply(object[["submodels"]], function(x){x[[1]][["model"]][["global"]]})))
  } else {
    n_global_i <- k_endogen_i * 0
  }
  
  # Deterministic variables
  n_i <- unlist(lapply(object[["submodels"]],
                       function(x){x[[1]][["model"]][["n"]]}))
  n <- unique(n_i)
  if (length(n) > 1) {
    stop("Number of deterministic terms differs across sub-models.\n",
         "Feel free to send a feature request.")
  }
  if (n > 0) {
    det_names_i <- lapply(object[["submodels"]], function(x){x[[1]][["model"]][["deterministic"]]})
    det_names <- unique(unlist(det_names_i))
  } else {
    det_names <- NULL
  }
  
  # Number of draws
  draws_i <- unlist(lapply(object[["submodels"]],
                           function(x){nrow(x[[1]][["posterior"]][["u_sigma_inv"]][["coeffs"]])}))
  draws <- unique(draws_i)
  if (length(draws) > 1) {
    stop("Number of posterior draws differs across sub-models.\n",
         "Feel free to send a feature request.")
  }
  
  # Number of observations
  tt_i <- unlist(lapply(object[["submodels"]],
                        function(x){nrow(x[[1]][["data"]][["train"]][["y"]])}))
  tt <- unique(tt_i)
  if (length(tt) > 1) {
    stop("Number of available observations differs across sub-models.\n",
         "Feel free to send a feature request.")
  }
  
  # Last observations
  end_i <- unlist(lapply(object[["submodels"]],
                         function(x){stats::tsp(x[[1]][["data"]][["train"]][["y"]])[2]}))
  end <- unique(end_i)
  if (length(end) > 1) {
    stop("Last date of last period differs across sub-models.\n",
         "Feel free to send a feature request.")
  }
  
  # Create the bvarmodel object ----
  
  # Model
  model <- NULL
  model[["type"]] <- "GVAR"
  model[["k"]] <- k
  model[["p"]] <- p
  model[["m"]] <- m
  model[["s"]] <- s
  model[["n"]] <- n
  model[["varsel"]] <- "none"
  endogen <- unlist(endogen_i)
  names(endogen) <- NULL
  model[["endogen"]] <- endogen
  if (global) {
    model[["exogen"]] <- global_i
  }
  model[["deterministic"]] <- det_names
  model[["error"]] <- "Not available"
  model[["tvp"]] <- FALSE
  model[["structural"]] <- TRUE
  model[["iterations"]] <- as.integer(draws)
  model[["burnin"]] <- 0L
  
  # Data
  endogen_data <- object[["global"]][["endogen"]][, unlist(pos_i)]
  y_tsp <- stats::tsp(endogen_data)
  data <- stats::embed(endogen_data, p + 1)
  
  y <- stats::ts(data[, 1:k], end = y_tsp[2], frequency = y_tsp[3])
  dimnames(y)[[2]] <- endogen
  
  x <- stats::ts(data[, -(1:k)], end = y_tsp[2], frequency = y_tsp[3])
  x_names <- NULL
  for (i in 1:p) {
    x_names <- append(x_names, paste0(endogen, ".l", i))
  }
  
  if (global) {
    tsp_global <- stats::tsp(object[["global"]][["exogen"]])
    global_data <- stats::embed(object[["global"]][["exogen"]][, global_i], s + 1)
    global_data <- stats::ts(global_data, end = tsp_global[2], frequency = tsp_global[3])
    x <- stats::na.omit(cbind(x, global_data))
    for (i in 1:(s + 1)) {
      x_names <- append(x_names, paste0(global_i, ".l", i - 1)) 
    }
  } else {
    global_data <- NULL
  }
  
  if (n > 0) {
    if (n == 1) {
      if (all(det_names_i == "const")) {
        x <- cbind(x, 1)
        x_names <- append(x_names, "const")
        det_data <- stats::ts(matrix(rep(1, nrow(x))), end = y_tsp[2], frequency = y_tsp[3],
                              class = c("mts", "ts", "matrix"))
      } else {
        stop("Revise deterministic.")
      }
    } else {
      stop("Revise deterministic.")
    }
  } else {
    det_data <- NULL
  }
  dimnames(x)[[2]] <- x_names
  
  data <- list("original" = list("endogen" = endogen_data,
                                 "exogen" = global_data,
                                 "deterministic" = det_data),
               "train" = list("y" = y,
                              "x" = x))
  
  # Create skeleton
  posterior <- NULL
  posterior[["a0"]] <- matrix(NA_real_, draws, k * k)
  posterior[["a"]] <- matrix(NA_real_, draws, k * k * p)
  
  if (global) {
    posterior[["b"]] <- matrix(NA, draws, k * m * (1 + s))
  } else {
    posterior[["b"]] <- NULL
  }
  
  if (n > 0) {
    posterior[["c"]] <- matrix(NA, draws, k * n)
  } else {
    posterior[["c"]] <- NULL
  }
  
  u_sigma_inv <- matrix(NA, draws, k * k)
  
  tvp_i <- unlist(lapply(object[["submodels"]], function(x){x[[1]][["model"]][["tvp"]]}))
  tvp <- any(tvp_i)
  if (tvp) {
    stop("TVP functionality not implemented yet.\n",
         "Feel free to send a feature request.")
  }
  
  sv_i <- unlist(lapply(object[["submodels"]], function(x){x[[1]][["model"]][["error"]] %in% c("sv", "sv+covar")}))
  sv <- any(sv_i)
  if (sv) {
    stop("SV functionality not implemented yet.\n",
         "Feel free to send a feature request.")
  }
  
  period <- which(stats::time(object[["global"]][["endogen"]]) == end)
  # if (is.null(period)) {
  #   period <- tt
  # } else {
  #   if (period > tt | period < 1) {
  #     stop("Implausible specification of argument 'period'.")
  #   }
  # }
  
  # Weight matrix ----
  # Get a list of the used weight matrices
  w_positions <- .get_weight_matrix_positions(object)
  w <- list()
  for (i in submodels) {
    pos <- (period - 1) * w_positions[[i]][["n_z_old"]] + w_positions[[i]][["model"]]
    w[[i]] <- object[["weights"]][[i]][pos,]
  }
  names(w) <- submodels
  
  cat(paste("Combining submodels to global model...\n"))
  pb <- utils::txtProgressBar(style = 3)
  for (draw in 1:draws) {
    
    # Put together A0 ----
    a0_temp <- matrix(NA_real_, k, k)
    for (i in submodels) {
      
      pos_a0_exogen <- n_endogen_i[i] + 1:(k_endogen_i[i] * k_exogen_i[i])
      
      # A0
      if (object[["submodels"]][[i]][[1]][["model"]][["structural"]]) {
        stop("implement structural")
        if (!is.null(object[[i]][["posterior"]][["a"]][["coeffs"]]) & tvp_i[i, "a0"]) {
          A0 <- matrix(object[[i]][["posterior"]][["a"]][["coeffs"]][[period]][draw, ], k_endogen_i[i]) 
        } else {
          A0 <- matrix(object[[i]][["posterior"]][["a"]][["coeffs"]][draw, ], k_endogen_i[i]) 
        }
      } else {
        A0 <- diag(1, k_endogen_i[i])
      }
      
      # Contemporary foreign
      if (tvp_i[i]) {
        stop("Implement TVP")
        A0_exogen <- matrix(object[[i]][["posterior"]][["a"]][["coeffs"]][[period]][draw, pos_a0_exogen], k_endogen_i[i])
      } else {
        A0_exogen <- matrix(object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_a0_exogen], k_endogen_i[i])
      }
      
      a0_temp[k_pos[i] + 1:k_endogen_i[i],] <- cbind(A0, -A0_exogen) %*% w[[i]]
    }
    if (!all(diag(a0_temp) == 1)) {
      stop("Not all diagonal elements of the global A0 matrix are 1.")
    }
    posterior[["a0"]][draw, ] <- matrix(a0_temp)
    rm(a0_temp)
    
    # Put together A ----
    for (j in 1:p) {
      
      # Create global matrix of [A_d, A_*] for lag j
      # Fill with NA to facilitate error spotting
      a_temp <- matrix(NA_real_, k, k)
      
      for (i in submodels) {
        
        # Create a sub-model matrix [A_d, A_*] and fill it with A_d and A_*
        temp_i <- matrix(0, k_endogen_i[i], k_endogen_i[i] + k_exogen_i[i])
        
        # Endogenous draws of lag j
        if (j <= p_endogen_i[i]) { # If j is larger than p_endogen_i, leave the country matrix 0
          pos_endogen <- (j - 1) * k_endogen_i[i] * k_endogen_i[i] + 1:(k_endogen_i[i] * k_endogen_i[i])
          if (tvp_i[i]) {
            stop("Implement TVP")
            temp_i[, 1:k_endogen_i[i]] <- object[[i]][["posterior"]][["domestic"]][[period]][draw, pos_endogen]
          } else {
            temp_i[, 1:k_endogen_i[i]] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_endogen] 
          }
        }
        
        # Foreign draws of lag j
        if (j <= p_exogen_i[i]) { # If j is larger than p_exogen_i, leave the country matrix 0
          pos_exogen <-  n_endogen_i[i] + j * (k_endogen_i[i] * k_exogen_i[i]) + 1:(k_endogen_i[i] * k_exogen_i[i])
          if (tvp_i[i]) {
            stop("Implement TVP")
            temp_i[, k_endogen_i[i] + 1:k_exogen_i[i]] <- object[[i]][["posterior"]][["a"]][["coeffs"]][[period]][draw, pos_exogen]  
          } else {
            temp_i[, k_endogen_i[i] + 1:k_exogen_i[i]] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_exogen] 
          }
        }
        
        a_temp[k_pos[i] + 1:k_endogen_i[i],] <- temp_i %*% w[[i]]
        rm(temp_i)
      }
      
      # Store
      posterior[["a"]][draw, (j - 1) * k^2 + 1:(k^2)] <- matrix(a_temp)
      rm(a_temp)
    }
    
    # Put together H ----
    if (global) {
      for (j in 1:(s + 1)) {
        h_temp <- matrix(0, k, m)
        for (i in submodels) {
          if (m_i[i] > 0) {
            if (j <= s_i[i] + 1) {
              pos_global <- n_endogen_i[i] + n_exogen_i[i] + (j - 1) * k_endogen_i[i] * m_i[i] + 1:(k_endogen_i[i] * m_i[i])
              if (tvp_i[i]) {
                stop("Implement TVP")
                #h_temp[k_pos[i] + 1:k_endogen_i[i],] <- object[[i]][["posterior"]][["a"]][["coeffs"]][[period]][draw, pos_global]
              } else {
                h_temp[k_pos[i] + 1:k_endogen_i[i],] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_global] 
              }
            }  
          }
        }
        
        # Store
        posterior[["b"]][draw, (j - 1) * k * m + 1:(k * m)] <- matrix(h_temp)
        rm(h_temp)
      }
    }
    
    #### Put together D ####
    if (n > 0) {
      d_temp <- matrix(0, k, n)
      for (i in submodels) {
        pos_det <- n_endogen_i[i] + n_exogen_i[i] + n_global_i[i] + 1:(k_endogen_i[i] * n_i[i])
        if (n_i[i] > 0) {
          if (tvp_i[i]) {
            stop("Implement TVP")
            #d_temp[k_pos[i] + 1:k_endogen_i[i],] <- object[[i]][["posterior"]][["a"]][["coeffs"]][[period]][draw, pos_det]
          } else {
            d_temp[k_pos[i] + 1:k_endogen_i[i],] <- object[["submodels"]][[i]][[1]][["posterior"]][["a"]][["coeffs"]][draw, pos_det] 
          }
        }
      }
      
      # Premultiply by A0_i and store
      posterior[["c"]][draw, ] <- matrix(d_temp)
      rm(d_temp)
    }
    
    # Put together Sigma ----
    sigma_temp <- matrix(0, k, k)
    for (i in submodels) {
      if (sv_i[i]) {
        stop("Implement SV")
        #sigma_temp[k_pos[i] + 1:k_endogen_i[i], k_pos[i] + 1:k_endogen_i[i]] <- object[[i]][["posterior"]][["sigma"]][["coeffs"]][[period]][draw, ] 
      } else {
        sigma_temp[k_pos[i] + 1:k_endogen_i[i], k_pos[i] + 1:k_endogen_i[i]] <- object[["submodels"]][[i]][[1]][["posterior"]][["u_sigma_inv"]][["coeffs"]][draw, ]
      }
    }
    u_sigma_inv[draw, ] <- sigma_temp
    rm(sigma_temp)
    
    utils::setTxtProgressBar(pb, value = draw / draws)
  }
  
  # Posterior
  coeffs <- coda::mcmc(cbind(posterior[["a"]], posterior[["b"]], posterior[["c"]], posterior[["a0"]]))
  u_sigma_inv <- coda::mcmc(u_sigma_inv)
  
  result <- list("model" = model,
                 "data" = data,
                 "posterior" = list("a" = list("coeffs" = coeffs),
                                    "u_sigma_inv" = list("coeffs" = u_sigma_inv)))
  
  class(result) <- list("bvarmodel", "list")
  return(result)
}