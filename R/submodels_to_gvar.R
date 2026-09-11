#' Generate Global Model from Sub-Models
#'
#' Combines the estimated sub-models of a global model into the solved global
#' model.
#'
#' @param object an object of class 'gvarmodel' containing exactly one
#' estimated model per sub-model. Use the \code{submodels} argument of
#' \code{\link{read_gvar_from_folder}} to reduce an object with multiple
#' candidate models per sub-model to one.
#' @param period integer of the period, whose weight matrices should be used to
#' solve the model. Defaults to the last period of the estimation sample. Only
#' relevant for time varying weights.
#'
#' @details
#' Each sub-model \eqn{i} is a VARX* model in its endogenous variables
#' \eqn{x_{it}} and its weakly exogenous variables \eqn{x^{*}_{it}}, which
#' are weighted averages of the endogenous variables of the remaining
#' sub-models. Written in terms of the vector \eqn{y_t} of the endogenous
#' variables of all sub-models, the contemporaneous terms of sub-model \eqn{i}
#' become \eqn{(A_{0i}, -A^{*}_{0i}) W_i y_t}, and stacking these blocks over
#' all sub-models gives the global matrix \eqn{G}. The lagged terms are
#' stacked the same way, which turns the sub-models into a single VAR in
#' \eqn{y_t}.
#'
#' The model is returned in its reduced form, that is, premultiplied by
#' \eqn{G^{-1}}, which is the form the methods of \pkg{bvartools} work on.
#' \eqn{G} is a dense matrix and can therefore not be passed on as the
#' structural form of a 'bvarmodel', which stores the free elements of a unit
#' lower triangular matrix. Its draws are returned as element \code{g}
#' instead.
#'
#' The error covariance matrix of the structural form is block diagonal, one
#' block per sub-model. The reduced form errors \eqn{G^{-1} u_t} are
#' correlated across sub-models accordingly.
#'
#' @return An object of class 'bvarmodel', which can be analysed with the
#' functionality of \pkg{bvartools}.
#'
#' @examples
#'
#' # Load data
#' data("gvar2019")
#' global_data <- gvar2019[["global_data"]]
#' submodel_data <- gvar2019[["submodel_data"]]
#'
#' # Limit number of sub-models
#' submodel_data <- select_list_elements(submodel_data, c("AT", "DE", "US"))
#'
#' # Set up the model
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = global_data)
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data,
#'                               period = 3)
#' object <- add_submodels(object,
#'                         p_endogen = 1,
#'                         exogen = c("y", "Dp", "eq", "r", "lr"),
#'                         p_exogen = 1,
#'                         global = "poil", s = 1,
#'                         error = "wishart",
#'                         iterations = 20, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' object <- align_model_obs(object)
#' object <- add_priors(object,
#'                      coef = list(v_i = 0),
#'                      sigma = list(df = 3, scale = 0.0001))
#' object <- add_initial_values(object)
#' object <- add_posterior_coefficients(object)
#'
#' # Solve the global model
#' gvar <- submodels_to_gvar(object)
#'
#' @export
submodels_to_gvar <- function(object, period = NULL) {
  
  index <- object[["global"]][["index"]]
  submodels <- unique(index[, "submodel"])

  # The sub-models are stacked as they are, so they have to be VARX models
  # already. A VECX sub-model has to be rewritten in levels first.
  classes <- unlist(lapply(object[["submodels"]],
                           function(x) class(x[[1]])))
  if (any(classes == "vecxsubmodel")) {
    stop("Argument 'object' contains sub-models in error correction form. ",
         "Use 'gvec_to_gvar' to obtain their level representation first.")
  }

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
    # match() the other way round would give, for every row of index_i, its
    # position among the model's variables -- not the row of index_i that
    # belongs to the j-th variable of the model, which is what is wanted here.
    pos_i[[i]] <- index_i[match(endogen_i[[i]], index_i[, "variable"]), "id"]
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

  # Without a lag nothing propagates between the units, so there is no dynamic
  # global model to solve.
  if (p == 0) {
    stop("None of the sub-models uses a lag of its endogenous or of its weakly ",
         "exogenous variables, so the global model has no dynamics.")
  }
  # Position of the first variable of each sub-model in the global vector.
  # Written without indexing so that a global model of a single sub-model,
  # where the shift would be an empty assignment, works as well.
  k_pos <- cumsum(c(0, k_endogen_i[-length(k_endogen_i)]))
  names(k_pos) <- names(k_endogen_i)
  
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
    s <- 0
    global_i <- NULL
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
  # The global model is returned in its reduced form (see below), whose errors
  # have a constant covariance matrix.
  model[["error"]] <- "wishart"
  model[["tvp"]] <- FALSE
  # The contemporaneous matrix G is inverted out below, so what is returned is
  # a reduced form model. G itself is a dense matrix, and the structural form
  # of 'bvartools' stores a unit lower triangular matrix of k * (k - 1) / 2
  # free elements, which cannot hold it.
  model[["structural"]] <- FALSE
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

  # The global model has to describe the sample the sub-models were estimated
  # on. Embedding the global series only drops the leading observations that
  # the lag structure of the global model needs, which can leave one more than
  # align_model_obs left the sub-models with -- for instance when the models
  # that were selected use fewer lags than the longest candidate did.
  sub_tsp <- stats::tsp(object[["submodels"]][[submodels[1]]][[1]][["data"]][["train"]][["y"]])
  if (sub_tsp[1] < stats::tsp(y)[1] - 1e-6) {
    stop("The sub-models were estimated from ", sub_tsp[1], " on, which the ",
         "lag structure of the global model leaves no room for.")
  }
  y <- stats::window(y, start = sub_tsp[1], end = sub_tsp[2])
  x <- stats::window(x, start = sub_tsp[1], end = sub_tsp[2])

  if (n > 0) {

    if (length(det_names) != n) {
      stop("The sub-models use ", n, " deterministic term(s) each, but ",
           length(det_names), " different ones between them: ",
           paste0(det_names, collapse = ", "), ".")
    }

    # From the global model, which is where the sub-models took them from as
    # well, so the global model uses the same series they were estimated with.
    det_data <- .submodel_deterministic(object, det_names,
                                        start = sub_tsp[1], end = sub_tsp[2])

    if (nrow(det_data) != nrow(x)) {
      stop("The deterministic terms span ", nrow(det_data), " observations, ",
           "the global model ", nrow(x), ".")
    }

    x <- cbind(x, det_data)
    x_names <- append(x_names, det_names)

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
  
  # Period of the weight matrices ----
  #
  # Defaults to the last period of the estimation sample. Compared with a
  # tolerance rather than for equality, since the times of a 'ts' object are
  # accumulated in floating point and a quarter of a year is not exact.
  times <- stats::time(object[["global"]][["endogen"]])

  if (is.null(period)) {
    period <- which.min(abs(times - end))
    if (abs(times[period] - end) > 1e-6) {
      stop("The last period of the sub-models is not part of the global data.")
    }
  } else {
    if (length(period) != 1 || period < 1 || period > length(times)) {
      stop("Argument 'period' must be a single integer between 1 and ",
           length(times), ".")
    }
  }

  # Constant weights are stored as a single matrix, which is then used for
  # every period.
  n_periods_w <- vapply(submodels, function(i) {
    n_vars <- length(unique(index[index[, "submodel"] == i, "variable"])) +
      length(unique(index[index[, "submodel"] != i, "variable"]))
    nrow(object[["weights"]][[i]]) / n_vars
  }, numeric(1))
  period_w <- if (all(n_periods_w == 1)) 1 else period

  # Weight matrix ----
  #
  # Rows in the order the coefficients of a sub-model are in, columns
  # restricted to the variables that are endogenous to the global model.
  keep <- unlist(pos_i)
  names(keep) <- NULL
  w <- lapply(submodels, .submodel_weight_matrix, object = object,
              period = period_w, keep = keep)
  names(w) <- submodels
  
  message("Combining sub-models to the global model ...")

  # Only interactively: a progress bar in the output of R CMD check or of a
  # batch job is noise.
  pb <- NULL
  if (interactive()) {
    pb <- utils::txtProgressBar(style = 3)
    on.exit(close(pb), add = TRUE)
  }
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
    # A sub-model contributes its own variables through an identity block of
    # its weight matrix, so the diagonal must come out as one. Compared with a
    # tolerance, since it is the result of a matrix product.
    if (any(abs(diag(a0_temp) - 1) > 1e-8)) {
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
    
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, value = draw / draws)
    }
  }
  
  # Reduced form ----
  #
  # Stacking the sub-models gives the structural form
  #   G y_t = A_1 y_{t-1} + ... + B x_t + C d_t + u_t,
  # in which G is the dense matrix of contemporaneous coefficients assembled
  # above. Every method that follows works on the reduced form
  #   y_t = G^-1 A_1 y_{t-1} + ... + G^-1 B x_t + G^-1 C d_t + G^-1 u_t,
  # so the model is solved for it here rather than leaving the inversion to
  # each of them.
  #
  # G is not triangular, so it cannot be passed on as the structural form of a
  # 'bvarmodel', which stores the k * (k - 1) / 2 free elements of a unit lower
  # triangular matrix. Its draws are returned separately, as element 'g'.
  for (draw in 1:draws) {

    g_i <- matrix(posterior[["a0"]][draw, ], k)
    g_inv <- tryCatch(solve(g_i), error = function(e) {
      stop("The global matrix of contemporaneous coefficients is singular in ",
           "draw ", draw, ", so the global model cannot be solved.")
    })

    for (j in seq_len(p)) {
      pos <- (j - 1) * k^2 + 1:(k^2)
      posterior[["a"]][draw, pos] <- g_inv %*% matrix(posterior[["a"]][draw, pos], k)
    }

    if (global) {
      posterior[["b"]][draw, ] <- g_inv %*% matrix(posterior[["b"]][draw, ], k)
    }

    if (n > 0) {
      posterior[["c"]][draw, ] <- g_inv %*% matrix(posterior[["c"]][draw, ], k)
    }

    # The reduced form error is G^-1 u_t, so its covariance is
    # G^-1 Sigma_u G^-1', and the precision that is stored is G' Sigma_u^-1 G.
    # Written in terms of the precision, so that the block diagonal Sigma_u
    # never has to be inverted.
    omega <- matrix(u_sigma_inv[draw, ], k)
    u_sigma_inv[draw, ] <- t(g_i) %*% omega %*% g_i
  }

  # Posterior
  coeffs <- coda::mcmc(cbind(posterior[["a"]], posterior[["b"]], posterior[["c"]]))
  u_sigma_inv <- coda::mcmc(u_sigma_inv)

  result <- list("model" = model,
                 "data" = data,
                 "posterior" = list("a" = list("coeffs" = coeffs),
                                    "u_sigma_inv" = list("coeffs" = u_sigma_inv)),
                 "g" = coda::mcmc(posterior[["a0"]]))

  class(result) <- list("bvarmodel", "list")
  return(result)
}