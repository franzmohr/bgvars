#' Generate a GVAR Model
#' 
#' Combines the country model results to a global VAR model and solves it.
#' 
#' @param object a list containing the results of the country estimates, usually, the
#' result of a call to \code{\link{draw_posterior}}.
#' @param period an integer of the time index for which the GVAR should be solved. Only used
#' when time varying weights or parameters are used.
#' @param thin an integer specifying the thinning factor for the MCMC output.
#' Defaults to 1 to obtain the full MCMC sequence.
#' 
#' @return An object of class \code{"bgvar"}.
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
combine_submodels <- function(object, period = NULL, thin = 1){
  
  if ("bgvecest" %in% class(object)) {
    stop("Please transform the country-specific VECX models to VARs before using this function.")
  }
  
  # Check if only one model per country
  names_object <- names(object)
  for (i in names_object) {
    if (sum(names_object == i) > 1) {
      stop("Multiple models for the same country detected. Please choose final models before using this function.")
    }
  }
  
  #### Solve GVAR model ####
  
  # Obtain number of draws per country model
  draws_i <- unlist(lapply(object, function(x){ifelse(class(x[["posteriors"]][["foreign"]]) == "list",
                                                      nrow(x[["posteriors"]][["foreign"]][[1]]),
                                                      nrow(x[["posteriors"]][["foreign"]]))}))
  
  # Check if the number of posterior draws is equal across country models
  draws <- unique(draws_i)
  if (length(draws) > 1) {
    stop("Number of posterior draws must be equal across all country models.")
  }
  
  tt_i <- unlist(lapply(object, function(x) {NROW(x[["data"]][["Y"]])}))
  tt <- unique(tt_i)
  if (length(tt) > 1) {
    stop("Number of observations is not equal across all country models.")
  }
  
  # Get variables for each country model
  n_countries <- length(object)
  country_names <- names(object)
  
  k_domestic_i <- unlist(lapply(object, function(x){length(x[["model"]][["domestic"]][["variables"]])}))
  k_foreign_i <- unlist(lapply(object, function(x){length(x[["model"]][["foreign"]][["variables"]])}))
  k <- sum(k_domestic_i) # Total number of endogenous variables
  k_pos <- cumsum(k_domestic_i) - k_domestic_i
  # Get lags for each country model
  p_domestic_i <- unlist(lapply(object, function(x){x[["model"]][["domestic"]][["lags"]]}))
  p_foreign_i <- unlist(lapply(object, function(x){x[["model"]][["foreign"]][["lags"]]}))
  p <- max(c(p_domestic_i, p_foreign_i)) # Lag of endogenous variables in the global model
  
  # Produce a vector of positions. Each position will be used during solving
  # the model while the others will be omitted.
  pos_thin <- seq(from = thin, to = draws, by = thin)
  draws <- length(pos_thin)
  
  # Create skeleton
  result <- NULL
  result[["a0"]] <- matrix(NA, k^2, draws)
  result[["a"]] <- matrix(NA, k^2 * p, draws)
  
  global_i <- unlist(lapply(object, function(x){!is.null(x[["data"]][["global"]])}))
  global <- any(global_i)
  if (global) {
    # Get global variables
    k_global_i <- unlist(lapply(object, function(x){length(x[["model"]][["global"]][["variables"]])}))
    k_global <- max(k_global_i)
    p_global_i <- unlist(lapply(object, function(x){x[["model"]][["global"]][["lags"]]}))
    s <- max(p_global_i) # Lag of global variables in the global model
    result[["b"]] <- matrix(0, k * k_global * (1 + s), draws)
  }
  
  # Get deterministic terms
  deter_i <- unlist(lapply(object, function(x){!is.null(x[["data"]][["deterministic"]])}))
  deter <- any(deter_i)
  if (deter) {
    k_deter_i <- unlist(lapply(object, function(x){ifelse(class(x[["data"]][["deterministic"]]) == "list",
                                                          ncol(x[["data"]][["deterministic"]][[1]]),
                                                          ncol(x[["data"]][["deterministic"]]))}))
    k_deter <- max(k_deter_i)
    result[["c"]] <- matrix(NA, k * k_deter, draws)
  }
  
  result[["sigma"]] <- matrix(NA, k * k, draws)
  
  tvp_i <- matrix(NA, n_countries, 6)
  dimnames(tvp_i) <- list(country_names,
                          c("domestic", "foreign", "global", "deterministic", "a0", "sigma"))
  for (i in country_names) {
    for (j in c("domestic", "foreign", "global", "deterministic", "a0", "sigma")) {
      if (!is.null(object[[i]][["posteriors"]][[j]])) {
        tvp_i[i, j] <- is.list(object[[i]][["posteriors"]][[j]]) 
      } 
    }
  }
  tvp <- any(tvp_i, na.rm = TRUE)
  
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
  
  cat(paste("Generating GVAR model...\n"))
  pb <- utils::txtProgressBar(style = 3)
  for (draw_i in 1:draws) {
    
    # Select draw that should be used considering thinning
    draw <- pos_thin[draw_i]
    
    #### Put together A0 ####
    a0_temp <- NULL
    for (i in country_names) {
      
      # Structural
      if (object[[i]][["model"]][["structural"]]) {
        if (!is.null(object[[i]][["posteriors"]][["a0"]]) & tvp_i[i, "a0"]) {
          A0 <- matrix(object[[i]][["posteriors"]][["a0"]][[period]][draw_i, ], k_domestic_i[i]) 
        } else {
          A0 <- matrix(object[[i]][["posteriors"]][["a0"]][draw_i, ], k_domestic_i[i]) 
        }
      } else {
        A0 <- diag(1, k_domestic_i[i])
      }
      
      # Contemporary foreign
      if (tvp_i[i, "foreign"]) {
        A0_for <- matrix(object[[i]][["posteriors"]][["foreign"]][[period]][draw_i, 1:(k_domestic_i[i] * k_foreign_i[i])], k_domestic_i[i])
      } else {
        A0_for <- matrix(object[[i]][["posteriors"]][["foreign"]][draw_i, 1:(k_domestic_i[i] * k_foreign_i[i])], k_domestic_i[i])
      }
      
      a0_temp <- rbind(a0_temp, cbind(A0, -A0_for) %*% W[[i]])
    }
    result[["a0"]][, draw_i] <- a0_temp
    
    #### Put together G ####
    for (j in 1:p) {
      # Create global matrix of [A_d, A_*] for lag j
      g_temp <- matrix(NA_real_, k, k)
      for (i in country_names) {
        # Create a country matrix [A_d, A_*] and fill it with A_D and A_*
        temp_i <- matrix(0, k_domestic_i[i], k_domestic_i[i] + k_foreign_i[i])
        
        # Domestic draws of lag j
        if (j <= p_domestic_i[i]) { # If j is larger than p_domestic_i, leave the country matrix 0
          if (tvp_i[i, "domestic"]) {
            temp_i[, 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["domestic"]][[period]][draw, (j - 1) * k_domestic_i[i]^2 + 1:k_domestic_i[i]^2]
          } else {
            temp_i[, 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["domestic"]][draw, (j - 1) * k_domestic_i[i]^2 + 1:k_domestic_i[i]^2] 
          }
        }
        # foreign draws of lag j
        if (tvp_i[i, "foreign"]) {
          if (j <= p_foreign_i[i]) { # If j is larger than p_foreign_i, leave the country matrix 0
            temp_i[, k_domestic_i[i] + 1:k_foreign_i[i]] <- object[[i]][["posteriors"]][["foreign"]][[period]][draw, j * k_domestic_i[i] * k_foreign_i[i]  + 1:(k_domestic_i[i] * k_foreign_i[i])]  
          } 
        } else {
          if (j <= p_foreign_i[i]) { # If j is larger than p_foreign_i, leave the country matrix 0
            temp_i[, k_domestic_i[i] + 1:k_foreign_i[i]] <- object[[i]][["posteriors"]][["foreign"]][draw, j * k_domestic_i[i] * k_foreign_i[i]  + 1:(k_domestic_i[i] * k_foreign_i[i])]  
          }  
        }
        
        g_temp[k_pos[i] + 1:k_domestic_i[i],] <- temp_i %*% W[[i]]
        rm(temp_i)
      }
      
      # Store
      result[["a"]][(j - 1) * k^2 + 1:k^2, draw_i] <- matrix(g_temp)
      rm(g_temp)
    }
    
    #### Put together H ####
    if (global) {
      for (j in 1:(s + 1)) {
        h_temp <- matrix(0, k, k_global)
        for (i in country_names) {
          if (j <= p_global_i[i] + 1) {
            if (tvp_i[i, "global"]) {
              h_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["global"]][[period]][draw, (j - 1) * k_domestic_i[i] * k_global_i[i]  + 1:(k_domestic_i[i] * k_global_i[i])]
            } else {
              h_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["global"]][draw, (j - 1) * k_domestic_i[i] * k_global_i[i]  + 1:(k_domestic_i[i] * k_global_i[i])] 
            }
          } 
        }
        
        # Store
        result[["b"]][(j - 1) * k * k_global + 1:(k * k_global), draw_i] <- h_temp
        rm(h_temp)
      }
    }
    
    #### Put together D ####
    if (deter) {
      d_temp <- matrix(0, k, k_deter)
      for (i in country_names) {
        if (deter_i[i]) {
          if (tvp_i[i, "deterministic"]) {
            d_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["deterministic"]][[period]][draw, ]
          } else {
            d_temp[k_pos[i] + 1:k_domestic_i[i],] <- object[[i]][["posteriors"]][["deterministic"]][draw, ] 
          }
        }
      }
      
      # Premultiply by A0_i and store
      result[["c"]][, draw_i] <- d_temp
      rm(d_temp)
    }
    
    #### Put together Sigma ####
    sigma_temp <- matrix(0, k, k)
    for (i in country_names) {
      if (tvp_i[i, "sigma"]) {
        sigma_temp[k_pos[i] + 1:k_domestic_i[i], k_pos[i] + 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["sigma"]][[period]][draw, ] 
      } else {
        sigma_temp[k_pos[i] + 1:k_domestic_i[i], k_pos[i] + 1:k_domestic_i[i]] <- object[[i]][["posteriors"]][["sigma"]][draw, ]  
      }
    }
    result[["sigma"]][, draw_i] <- sigma_temp
    rm(sigma_temp)
    
    utils::setTxtProgressBar(pb, value = draw_i / draws)
  }
  
  # Convert posterior draws to coda objects ----
  mc_start <- unique(unlist(lapply(object, function(x){ifelse(class(x[["posteriors"]][["foreign"]]) == "list",
                                                              attributes(x[["posteriors"]][["foreign"]][[1]])$mcpar[1],
                                                              attributes(x[["posteriors"]][["foreign"]])$mcpar[1])})))
  mc_end <- unique(unlist(lapply(object, function(x){ifelse(class(x[["posteriors"]][["foreign"]]) == "list",
                                                            attributes(x[["posteriors"]][["foreign"]][[1]])$mcpar[2],
                                                            attributes(x[["posteriors"]][["foreign"]])$mcpar[2])})))
  mc_thin <- unique(unlist(lapply(object, function(x){ifelse(class(x[["posteriors"]][["foreign"]]) == "list",
                                                             attributes(x[["posteriors"]][["foreign"]][[1]])$mcpar[3],
                                                             attributes(x[["posteriors"]][["foreign"]])$mcpar[3])})))
  
  mc_seq <- seq(from = mc_start, to = mc_end, by = mc_thin)[pos_thin]
  mc_start <- mc_seq[1]
  mc_thin <- mc_seq[2] - mc_seq[1]
  
  result[["a0"]] <- coda::mcmc(t(result[["a0"]]), start = mc_start, thin = mc_thin)
  result[["a"]] <- coda::mcmc(t(result[["a"]]), start = mc_start, thin = mc_thin)
  if (global) {
    result[["b"]] <- coda::mcmc(t(result[["b"]]), start = mc_start, thin = mc_thin)
  }
  if (deter) {
    result[["c"]] <- coda::mcmc(t(result[["c"]]), start = mc_start, thin = mc_thin)
  }
  result[["sigma"]] <- coda::mcmc(t(result[["sigma"]]), start = mc_start, thin = mc_thin)
  
  #### Create variable index ----
  result[["data"]][["endogen"]] <- NULL
  data_names <- NULL
  country_names <- NULL
  for (i in names(object)) {
    result[["data"]][["endogen"]] <- cbind(result[["data"]][["endogen"]], object[[i]][["data"]][["domestic"]])
    data_names <- c(data_names, dimnames(object[[i]][["data"]][["domestic"]])[[2]])
    country_names <- c(country_names, rep(i, length(dimnames(object[[i]][["data"]][["domestic"]])[[2]])))
  }
  dimnames(result[["data"]][["endogen"]])[[2]] <- data_names
  index <- data.frame("country" = country_names,
                      "variable" = data_names,
                      stringsAsFactors = FALSE)
  result[["index"]] <- index
  
  # Model specs
  result[["model"]] <- list()
  result[["model"]][["endogen"]] <- list(variables = data_names, lags = p)
  
  
  # Collect raw data ----
  ## Put together deterministic terms ----
  if (deter) {
    data_names <- NULL
    for (i in names(object)) {
      if (object[[i]]$model$type == "VEC" & !is.null(object[[i]]$data$deterministic)) {
        if (!is.null(object[[i]]$data$deterministic$unrestricted)) {
          result$data$deterministic <- cbind(result$data$deterministic, object[[i]]$data$deterministic$unrestricted)
          data_names <- c(data_names, dimnames(object[[i]]$data$deterministic$unrestricted)[[2]])
        }
        if (!is.null(object[[i]]$data$deterministic$restricted)) {
          result$data$deterministic <- cbind(result$data$deterministic, object[[i]]$data$deterministic$restricted) 
          data_names <- c(data_names, dimnames(object[[i]]$data$deterministic$restricted)[[2]])
        }
      } else {
        result$data$deterministic <- cbind(result$data$deterministic, object[[i]]$data$deterministic)
        data_names <- c(data_names, dimnames(object[[i]]$data$deterministic)[[2]])
      }
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
    
    # Update specs
    result[["model"]][["deterministic"]] <- list(variables = pos_det_name)
  }
  
  ## Put together global data ----
  if (global) {
    exogen <- NULL
    exogen_names <- NULL
    for (i in names(object)) {
      if (global_i[i]) {
        data_names_i <- dimnames(object[[i]][["data"]][["global"]])[[2]]
        pos <- which(!data_names_i %in% exogen_names)
        if (length(pos) > 0) {
          exogen <- cbind(exogen, object[[i]][["data"]][["global"]][, pos])
          exogen_names <- c(exogen_names, data_names_i[pos])        
        }
      }
    }
    temp_tsp <- stats::tsp(exogen)
    result[["data"]][["global"]] <- stats::as.ts(as.matrix(exogen))
    dimnames(result[["data"]][["global"]])[[2]] <- exogen_names
    stats::tsp(result[["data"]][["global"]]) <- temp_tsp
  
    # Update specs
    result[["model"]][["global"]] <- list(variables = exogen_names, lags = s)  
  }
  
  class(result) <- list("bgvar", "list")
  return(result)
}