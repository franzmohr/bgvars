#' Generate a GVAR Model
#' 
#' Combines the submodels of a global VAR model and solves it.
#' 
#' @param object a list containing the results of the submodel estimates, usually, the
#' result of a call to \code{\link{draw_posterior}}.
#' @param period an integer of the time index for which the GVAR should be solved. Only used
#' when time varying weights or parameters are used.
#' 
#' @return An object of class 'bgvar'.
#' 
#' @examples 
#' # Load data
#' data("gvar2019")
#' 
#' # Create regions
#' temp <- create_regions(country_data = gvar2019$country_data,
#'              weight_data = gvar2019$weight_data,
#'              region_weights = gvar2019$region_weights,
#'              regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")),
#'              period = 3)
#' 
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' global_data = gvar2019$global_data
#' 
#' # Difference series to make them stationary
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"), multi = 100)
#' global_data <- diff_variables(global_data, multi = 100)
#' 
#' # Create time varying weights
#' weight_data <- create_weights(weight_data, period = 3, country_data = country_data)
#' 
#' # Generate specifications
#' model_specs <- create_specifications(
#'                  country_data = country_data,
#'                  global_data = global_data,
#'                  countries = c("US", "JP", "CA", "NO", "GB", "EA"), 
#'                  domestic = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  foreign = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  global = list(variables = c("poil"), lags = 1),
#'                  deterministic = list(const = TRUE, trend = FALSE, seasonal = FALSE),
#'                  iterations = 10,
#'                  burnin = 10)
#' # Note that the number of iterations and burnin draws should be much higher!
#'                                      
#' # Overwrite country-specific specifications
#' model_specs[["US"]][["domestic"]][["variables"]] <- c("y", "Dp", "r")
#' model_specs[["US"]][["foreign"]][["variables"]] <- c("y", "Dp")
#' 
#' # Create estimation objects
#' country_models <- create_models(country_data = country_data,
#'                                 weight_data = weight_data,
#'                                 global_data = global_data,
#'                                 model_specs = model_specs)
#' 
#' # Add priors
#' models_with_priors <- add_priors(country_models,
#'                                  coef = list(v_i = 1 / 9, v_i_det = 1 / 10),
#'                                  sigma = list(df = 3, scale = .0001))
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(models_with_priors)
#' 
#' # Solve GVAR
#' gvar <- combine_submodels(object)
#' 
#' @export
combine_submodels <- function(object, period = NULL){
  
  # Check if only one model per country
  names_object <- names(object)
  if (any(table(names_object) > 1)) {
    stop("Argument 'object' contains more than one model for an entity or country.")
  }
  
  if ("vecxsubmodelest" %in% unlist(lapply(object, class))) {
    message("VECX submodels detected. Transforming them into their VARX representations...")
    for (i in 1:length(object)) {
      if ("vecxsubmodelest" %in% class(object[[i]])) {
        object[[i]] <- bvec_to_bvar(object[[i]])
      }
    }
  }
  
  #### Solve GVAR model ####
  
  # Obtain number of draws per country model
  # Use sigma variables, since those have to be always there
  draws_i <- unlist(lapply(object, function(x){ifelse(class(x[["posteriors"]][["sigma"]][["coeffs"]]) == "list",
                                                      nrow(x[["posteriors"]][["sigma"]][["coeffs"]][[1]]),
                                                      nrow(x[["posteriors"]][["sigma"]][["coeffs"]]))}))
  
  # Check if the number of posterior draws is equal across country models
  draws <- unique(draws_i)
  if (length(draws) > 1) {
    stop("Number of posterior draws must be equal across all country models.")
  }
  
  tt_i <- unlist(lapply(object, function(x) {NROW(x[["data"]][["y"]])}))
  tt <- unique(tt_i)
  if (length(tt) > 1) {
    stop("Number of observations is not equal across all country models.")
  }
  
  # Get variables for each country model
  n_countries <- length(object)
  country_names <- names(object)
  
  k_domestic_i <- unlist(lapply(object, function(x){x[["model"]][["k_domestic"]]}))
  k_foreign_i <- unlist(lapply(object, function(x){x[["model"]][["k_foreign"]]}))
  m_i <- unlist(lapply(object, function(x){x[["model"]][["m"]]}))
  m <- max(m_i)
  global <- m > 0
  n_i <- unlist(lapply(object, function(x){x[["model"]][["n"]]}))
  n <- max(n_i)
  k <- sum(k_domestic_i) # Total number of endogenous variables
  k_pos <- cumsum(k_domestic_i) - k_domestic_i
  # Get lags for each country model
  p_domestic_i <- unlist(lapply(object, function(x){x[["model"]][["p_domestic"]]}))
  p_foreign_i <- unlist(lapply(object, function(x){x[["model"]][["p_foreign"]]}))
  s_i <- unlist(lapply(object, function(x){x[["model"]][["s"]]}))
  p <- max(c(p_domestic_i, p_foreign_i)) # Lag of endogenous variables in the global model
  s <- max(s_i)
  
  n_domestic_i <- k_domestic_i * k_domestic_i * p_domestic_i
  n_foreign_i <- k_domestic_i * k_foreign_i * (p_foreign_i + 1)
  n_global_i <- k_domestic_i * m_i * (s_i + 1)
  n_c_i <- k_domestic_i * n_i
  
  # Create skeleton
  result <- NULL
  result[["a0"]] <- matrix(NA, draws, k * k)
  result[["a"]] <- matrix(NA, draws, k * k * p)
  
  
  if (global) {
    s <- max(s_i) # Lag of global variables in the global model
    result[["b"]] <- matrix(NA, draws, k * m * (1 + s))
  }
  
  # Get deterministic terms
  if (length(unique(n_i)) > 1) {
    stop("Heterogeneous specification of deterministic terms not supported yet.")
  }
  if (n > 0) {
    result[["c"]] <- matrix(NA, draws, k * n)
  }
  
  result[["sigma"]] <- matrix(NA, draws, k * k)
  
  tvp_i <- unlist(lapply(object, function(x){x[["model"]][["tvp"]]}))
  tvp <- any(tvp_i)
  
  sv_i <- unlist(lapply(object, function(x){x[["model"]][["error"]] %in% c("sv", "sv+covar")}))
  
  if (is.null(period)) {
    period <- tt
  } else {
    if (period > tt | period < 1) {
      stop("Implausible specification of argument 'period'.")
    }
  }
  
  #### Weight matrix ####
  # Get a list of the used weight matrices
  W <- NULL
  for (i in names(object)) {
    temp <- object[[i]]$data$weights
    if (!"matrix" %in% class(temp)) { # If weights are time varying...
      temp <- temp[,, period]
    }
    W <- c(W, list(temp))
  }
  names(W) <- names(object)
  rm(temp)
  
  cat(paste("Combining submodels to global model...\n"))
  pb <- utils::txtProgressBar(style = 3)
  for (draw in 1:draws) {
    
    #### Put together A0 ####
    a0_temp <- matrix(NA, k, k)
    for (i in country_names) {
      
      pos_a0_foreign <- n_domestic_i[i] + 1:(k_domestic_i[i] * k_foreign_i[i])
      # Structural
      if (object[[i]][["model"]][["structural"]]) {
        stop("implement structural")
        if (!is.null(object[[i]][["posteriors"]][["a"]][["coeffs"]]) & tvp_i[i, "a0"]) {
          A0 <- matrix(object[[i]][["posteriors"]][["a"]][["coeffs"]][[period]][draw, ], k_domestic_i[i]) 
        } else {
          A0 <- matrix(object[[i]][["posteriors"]][["a"]][["coeffs"]][draw, ], k_domestic_i[i]) 
        }
      } else {
        A0 <- diag(1, k_domestic_i[i])
      }
      
      # Contemporary foreign
      if (tvp_i[i]) {
        stop("Implement TVP")
        A0_for <- matrix(object[[i]][["posteriors"]][["a"]][["coeffs"]][[period]][draw, pos_a0_foreign], k_domestic_i[i])
      } else {
        A0_for <- matrix(object[[i]][["posteriors"]][["a"]][["coeffs"]][draw, pos_a0_foreign], k_domestic_i[i])
      }
      
      a0_temp[k_pos[i] + 1:k_domestic_i[i],] <- cbind(A0, -A0_for) %*% W[[i]]
    }
    result[["a0"]][draw, ] <- matrix(a0_temp)
    rm(a0_temp)
    
    #### Put together G ####
    for (j in 1:p) {
      # Create global matrix of [A_d, A_*] for lag j
      g_temp <- matrix(NA_real_, k, k)
      
      for (i in country_names) {
        
        # Create a country matrix [A_d, A_*] and fill it with A_D and A_*
        temp_i <- matrix(0, k_domestic_i[i], k_domestic_i[i] + k_foreign_i[i])
        
        # Domestic draws of lag j
        if (j <= p_domestic_i[i]) { # If j is larger than p_domestic_i, leave the country matrix 0
          pos_dom <- (j - 1) * k_domestic_i[i] * k_domestic_i[i] + 1:(k_domestic_i[i] * k_domestic_i[i])
          if (tvp_i[i]) {
            stop("Implement TVP")
            temp_i[, 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["domestic"]][[period]][draw, pos_dom]
          } else {
            temp_i[, 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["a"]][["coeffs"]][draw, pos_dom] 
          }
        }
        
        # Foreign draws of lag j
        if (j <= p_foreign_i[i]) { # If j is larger than p_foreign_i, leave the country matrix 0
          pos_for <-  n_domestic_i[i] + j * (k_domestic_i[i] * k_foreign_i[i]) + 1:(k_domestic_i[i] * k_foreign_i[i])
          if (tvp_i[i]) {
            stop("Implement TVP")
            temp_i[, k_domestic_i[i] + 1:k_foreign_i[i]] <- object[[i]][["posteriors"]][["a"]][["coeffs"]][[period]][draw, pos_for]  
          } else {
            temp_i[, k_domestic_i[i] + 1:k_foreign_i[i]] <- object[[i]][["posteriors"]][["a"]][["coeffs"]][draw, pos_for]  
          }
        }
        
        g_temp[k_pos[i] + 1:k_domestic_i[i],] <- temp_i %*% W[[i]]
        rm(temp_i)
      }
      
      # Store
      result[["a"]][draw, (j - 1) * k^2 + 1:(k^2)] <- matrix(g_temp)
      rm(g_temp)
    }
    
    #### Put together H ####
    if (global) {
      for (j in 1:(s + 1)) {
        h_temp <- matrix(0, k, m)
        for (i in country_names) {
          if (m_i[i] > 0) {
            if (j <= s_i[i] + 1) {
              pos_global <- n_domestic_i[i] + n_foreign_i[i] + (j - 1) * k_domestic_i[i] * m_i[i] + 1:(k_domestic_i[i] * m_i[i])
              if (tvp_i[i]) {
                h_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["a"]][["coeffs"]][[period]][draw, pos_global]
              } else {
                h_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["a"]][["coeffs"]][draw, pos_global] 
              }
            }  
          }
        }
        
        # Store
        result[["b"]][draw, (j - 1) * k * m + 1:(k * m)] <- matrix(h_temp)
        rm(h_temp)
      }
    }
    
    #### Put together D ####
    if (n > 0) {
      d_temp <- matrix(0, k, n)
      for (i in country_names) {
        pos_det <- n_domestic_i[i] + n_foreign_i[i] + n_global_i[i] + 1:(k_domestic_i[i] * n_i[i])
        if (n_i[i] > 0) {
          if (tvp_i[i]) {
            d_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["a"]][["coeffs"]][[period]][draw, pos_det]
          } else {
            d_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["a"]][["coeffs"]][draw, pos_det] 
          }
        }
      }
      
      # Premultiply by A0_i and store
      result[["c"]][draw, ] <- matrix(d_temp)
      rm(d_temp)
    }
    
    #### Put together Sigma ####
    sigma_temp <- matrix(0, k, k)
    for (i in country_names) {
      if (sv_i[i]) {
        sigma_temp[k_pos[i] + 1:k_domestic_i[i], k_pos[i] + 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["sigma"]][["coeffs"]][[period]][draw, ] 
      } else {
        sigma_temp[k_pos[i] + 1:k_domestic_i[i], k_pos[i] + 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["sigma"]][["coeffs"]][draw, ]  
      }
    }
    result[["sigma"]][draw, ] <- sigma_temp
    rm(sigma_temp)
    
    utils::setTxtProgressBar(pb, value = draw / draws)
  }
  
  # Convert posterior draws to coda objects ----
  mc_start <- unique(unlist(lapply(object, function(x){ifelse("list" %in% class(x[["posteriors"]][["a"]][["coeffs"]]),
                                                              attributes(x[["posteriors"]][["a"]][["coeffs"]][[1]])$mcpar[1],
                                                              attributes(x[["posteriors"]][["a"]][["coeffs"]])$mcpar[1])})))
  # mc_end <- unique(unlist(lapply(object, function(x){ifelse("list" %in% class(x[["posteriors"]][["a"]][["coeffs"]]),
  #                                                           attributes(x[["posteriors"]][["a"]][["coeffs"]][[1]])$mcpar[2],
  #                                                           attributes(x[["posteriors"]][["a"]][["coeffs"]])$mcpar[2])})))
  mc_thin <- unique(unlist(lapply(object, function(x){ifelse("list" %in% class(x[["posteriors"]][["a"]][["coeffs"]]),
                                                             attributes(x[["posteriors"]][["a"]][["coeffs"]][[1]])$mcpar[3],
                                                             attributes(x[["posteriors"]][["a"]][["coeffs"]])$mcpar[3])})))
  
  result[["a0"]] <- coda::mcmc(result[["a0"]], start = mc_start, thin = mc_thin)
  result[["a"]] <- coda::mcmc(result[["a"]], start = mc_start, thin = mc_thin)
  if (global) {
    result[["b"]] <- coda::mcmc(result[["b"]], start = mc_start, thin = mc_thin)
  }
  if (n > 0) {
    result[["c"]] <- coda::mcmc(result[["c"]], start = mc_start, thin = mc_thin)
  }
  result[["sigma"]] <- coda::mcmc(result[["sigma"]], start = mc_start, thin = mc_thin)
  
  # Collect raw data ----
  #### Create variable index ----
  result[["data"]][["y"]] <- NULL
  index <- NULL
  for (i in country_names) {
    result[["data"]][["y"]] <- cbind(result[["data"]][["y"]], object[[i]][["data"]][["domestic"]])
    data_names <- dimnames(object[[i]][["data"]][["domestic"]])[[2]]
    index <- rbind(index, data.frame("country" = i, "variable" = data_names, stringsAsFactors = FALSE))
  }
  dimnames(result[["data"]][["y"]])[[2]] <- index[, "variable"]
  
  ## Put together global data ----
  if (global) {
    exogen <- NULL
    exogen_names <- NULL
    for (i in country_names) {
      if (m_i[i]) {
        data_names_i <- dimnames(object[[i]][["data"]][["global"]])[[2]]
        pos <- which(!data_names_i %in% exogen_names)
        if (length(pos) > 0) {
          exogen <- cbind(exogen, object[[i]][["data"]][["global"]][, pos])
          exogen_names <- c(exogen_names, data_names_i[pos])        
        }
      }
    }
    temp_tsp <- stats::tsp(exogen)
    result[["data"]][["x"]] <- stats::as.ts(as.matrix(exogen))
    dimnames(result[["data"]][["x"]])[[2]] <- exogen_names
    stats::tsp(result[["data"]][["x"]]) <- temp_tsp
  }
  
  ## Put together deterministic terms ----
  if (n > 0) {
    data_names <- NULL
    for (i in country_names) {
      # if (object[[i]]$model$type == "VEC" & !is.null(object[[i]]$data$deterministic)) {
      #   if (!is.null(object[[i]]$data$deterministic$unrestricted)) {
      #     result$data$deterministic <- cbind(result$data$deterministic, object[[i]]$data$deterministic$unrestricted)
      #     data_names <- c(data_names, dimnames(object[[i]]$data$deterministic$unrestricted)[[2]])
      #   }
      #   if (!is.null(object[[i]]$data$deterministic$restricted)) {
      #     result$data$deterministic <- cbind(result$data$deterministic, object[[i]]$data$deterministic$restricted) 
      #     data_names <- c(data_names, dimnames(object[[i]]$data$deterministic$restricted)[[2]])
      #   }
      # } else {
      result$data$deterministic <- cbind(result$data$deterministic, object[[i]]$data$deterministic)
      data_names <- c(data_names, dimnames(object[[i]]$data$deterministic)[[2]])
      #}
    }
    dimnames(result$data$deterministic)[[2]] <- data_names
    pos_det_name <- NULL
    pos_det <- NULL
    for (i in 1:length(data_names)) {
      if (!data_names[i] %in% pos_det_name) {
        pos_det_name <- c(pos_det_name, data_names[i])
        pos_det <- c(pos_det, i)
      }
    }
    temp_tsp <- stats::tsp(result$data$deterministic)
    result$data$deterministic <- stats::as.ts(as.matrix(result$data$deterministic[, pos_det]))
    dimnames(result$data$deterministic)[[2]] <- pos_det_name
    stats::tsp(result$data$deterministic) <- temp_tsp
  }
  
  # Model specs
  result[["model"]] <- list()
  result[["model"]][["k"]] <- k
  result[["model"]][["p"]] <- p
  result[["model"]][["m"]] <- m
  result[["model"]][["s"]] <- s
  result[["model"]][["n"]] <- n
  result[["model"]][["index"]] <- index
  
  class(result) <- append("bgvar", class(result))
  return(result)
}
