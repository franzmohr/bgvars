#' Solve a GVAR Model
#' 
#' Combines the country model results to a global VAR model and solves it.
#' 
#' @param object a list containing the results of the country estimates, usually, the
#' result of a call to \code{\link{estimate_gvar}}.
#' @param ic a character specifying the information criterion used for model selection.
#' Available options are "AIC", "BIC" (default) and "HQ".
#' @param t an integer of the time index for which the GVAR should be solved. Only used
#' when time varying weights or parameters are used.
#' @param thin an integer specifying the thinning factor for the MCMC output.
#' Defaults to 1 to obtain the full MCMC sequence.
#' @param select a character specifying how the best country model is selected.
#' If \code{"order"} (default), the country model with the overall minimum value of the
#' specified information criterion per country is selected. If \code{"rank"}, the
#' function selects the model, after which the selected information crition increases
#' for the first time.
#' 
#' @details If multiple country models were estimated with \code{\link{estimate_gvar}},
#' the function will choose the model, which maximises the specified information
#' criterion specified by \code{ic}.
#' 
#' @return An object of class \code{"bgvar"}.
#' 
#' @examples 
#' data("gvar2016")
#' 
#' country_data <- gvar2016$country_data
#' global_data <- gvar2016$global_data
#' region_weights <- gvar2016$region_weights
#' weight_data <- gvar2016$weight_data
#' 
#' # Take first difference of all variables y and Dp
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"))
#' 
#' # Generate EA area region with 2 year rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 2,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 2 year rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "r"),
#'                                      countries = c("EA", "US", "JP", "CA", "MX", "GB"),
#'                                      p_domestic = 1,
#'                                      p_foreign = 1,
#'                                      type = "VAR")
#' 
#' # Create estimable objects
#' object <- create_models(country_data = country_data,
#'                         gvar_weights = gvar_weights,
#'                         model_specs = model_specs)
#' 
#' # Add priors
#' object <- add_priors(object)
#' 
#' # Estimate GVAR model
#' gvar_est <- estimate_gvar(object, iterations = 100, burnin = 10, thin = 2)
#' # Note that the number of iterations and burnin should be much higher, e.g., 30000.
#' 
#' # Solve GVAR
#' gvar_solved <- solve_gvar(gvar_est)
#' 
#' @export
solve_gvar <- function(object, ic = "BIC", t = NULL, thin = 1, select = "order"){
  #### Model selection ####
  
  criteria <- teststats(object, ic = ic, show_min = TRUE)
  data <- select_model(object, ic = ic, select = select)
  
  #### Solve GVAR model ####
  draws_i <- unlist(lapply(data, function(x){ifelse(class(x$coefficients$a_domestic) == "list",
                                                    nrow(x$coefficients$a_domestic[[1]]),
                                                    nrow(x$coefficients$a_domestic))}))
  
  draws <- unique(draws_i)
  if (length(draws) > 1) {
    stop("Number of posterior draws must be equal across all country models.")
  }
  
  k_domestic_i <- unlist(lapply(data, function(x){length(x$model$domestic$variables)}))
  k_foreign_i <- unlist(lapply(data, function(x){length(x$model$foreign$variables)}))
  k <- sum(k_domestic_i)
  k_pos <- cumsum(k_domestic_i) - k_domestic_i
  p_domestic_i <- unlist(lapply(data, function(x){x$model$domestic$lags}))
  p_foreign_i <- unlist(lapply(data, function(x){x$model$foreign$lags}))
  p <- max(c(p_domestic_i, p_foreign_i))
  
  pos_thin <- seq(from = thin, to = draws, by = thin)
  draws <- length(pos_thin)
  
  result <- NULL
  result$g0_i <- matrix(NA, k^2, draws)
  result$g <- matrix(NA, k^2 * p, draws)
  
  global_i <- unlist(lapply(data, function(x){!is.null(x$data$x_global)}))
  global <- any(global_i)
  if (global) {
    k_global_i <- unlist(lapply(data, function(x){length(x$model$global$variables)}))
    k_global <- max(k_global_i)
    p_global_i <- unlist(lapply(data, function(x){x$model$global$lags}))
    s <- max(p_global_i)
    result$h <- matrix(0, k * k_global * (1 + s), draws)
  }
  
  deter_i <- unlist(lapply(data, function(x){!is.null(x$data$deterministics)}))
  deter <- any(deter_i)
  if (deter) {
    k_deter_i <- unlist(lapply(data, function(x){ifelse(class(x$coefficients$c) == "list",
                                                        ncol(x$coefficients$c[[1]]),
                                                        ncol(x$coefficients$c))})) / k_domestic_i
    k_deter <- max(k_deter_i)
    result$d <- matrix(NA, k * k_deter, draws)
  }
  
  result$sigma <- matrix(NA, k * k, draws)
  
  #### Weight matrix ####
  W <- NULL
  for (i in names(data)) {
    temp <- data[[i]]$data$weights
    if (class(temp) == "array") {
      if (is.null(t)) {
        temp <- temp[,, dim(temp)[3]]
      } else {
        temp <- temp[,, t]
      }
    }
    W <- c(W, list(temp))
  }
  names(W) <- names(data)
  rm(temp)
  
  cat(paste("Solving GVAR model.\n"))
  for (draw_i in 1:draws) {
    draw <- pos_thin[draw_i]
    
    #### Put together A0 ####
    a0_temp <- NULL
    for (i in names(data)) {
        temp <- cbind(diag(1, k_domestic_i[i]),
                      matrix(-data[[i]]$coefficients$a_foreign[draw, 1:(k_domestic_i[i] * k_foreign_i[i])], k_domestic_i[i])) 

      temp <- temp %*% W[[i]]
      
      a0_temp <- rbind(a0_temp, temp)
      rm(temp)
    }
    a0_i_temp <- solve(a0_temp)
    
    #### G0_i ####
    result$g0_i[, draw_i] <- a0_i_temp
    
    #### Put together G ####
    for (j in 1:p) {
      g_temp <- matrix(NA, k, k)
      for (i in names(data)) {
        temp_i <- matrix(0, k_domestic_i[i], k_domestic_i[i] + k_foreign_i[i])
    
          if (j <= p_domestic_i[i]) {
            temp_i[, 1:k_domestic_i[i]] <- data[[i]]$coefficients$a_domestic[draw, (j - 1) * k_domestic_i[i]^2 + 1:k_domestic_i[i]^2]
          }
          if (j <= p_foreign_i[i]) {
            temp_i[, k_domestic_i[i] + 1:k_foreign_i[i]] <- data[[i]]$coefficients$a_foreign[draw, j * k_domestic_i[i] * k_foreign_i[i]  + 1:(k_domestic_i[i] * k_foreign_i[i])]  
          } 
        
        g_temp[k_pos[i] + 1:k_domestic_i[i],] <- temp_i %*% W[[i]]
        rm(temp_i)
      }
      result$g[(j - 1) * k^2 + 1:k^2, draw_i] <- a0_i_temp %*% g_temp
      rm(g_temp)
    }
    
    #### Put together H ####
    if (global) {
      for (j in 1:(s + 1)) {
        h_temp <- matrix(0, k, k_global)
        for (i in names(data)) {
          if (j <= p_global_i[i] + 1) {
            h_temp[k_pos[i] + 1:k_domestic_i[i],] <- data[[i]]$coefficients$a_global[draw, (j - 1) * k_domestic_i[i] * k_global_i[i]  + 1:(k_domestic_i[i] * k_global_i[i])]
          } 
        }
        result$h[(j - 1) * k * k_global + 1:(k * k_global), draw_i] <- a0_i_temp %*% h_temp
        rm(h_temp)
      }
    }
    
    #### Put together D ####
    if (deter) {
      d_temp <- matrix(0, k, k_deter)
      for (i in names(data)) {
          if (deter_i[i]) {
            d_temp[k_pos[i] + 1:k_domestic_i[i],] <- data[[i]]$coefficients$c[draw, ]
          }
      }
      result$d[, draw_i] <- a0_i_temp %*% d_temp
      rm(d_temp)
    }
    
    #### Put together Sigma ####
    sigma_temp <- matrix(0, k, k)
    for (i in names(data)) {
        sigma_temp[k_pos[i] + 1:k_domestic_i[i], k_pos[i] + 1:k_domestic_i[i]] <- data[[i]]$coefficients$sigma[draw, ] 
    }
    result$sigma[, draw_i] <- sigma_temp
    rm(sigma_temp)
  }
  
  mc_start <- unique(unlist(lapply(data, function(x){ifelse(class(x$coefficients$a_domestic) == "list",
                                                            attributes(x$coefficients$a_domestic[[1]])$mcpar[1],
                                                            attributes(x$coefficients$a_domestic)$mcpar[1])})))
  mc_end <- unique(unlist(lapply(data, function(x){ifelse(class(x$coefficients$a_domestic) == "list",
                                                          attributes(x$coefficients$a_domestic[[1]])$mcpar[2],
                                                          attributes(x$coefficients$a_domestic)$mcpar[2])})))
  mc_thin <- unique(unlist(lapply(data, function(x){ifelse(class(x$coefficients$a_domestic) == "list",
                                                           attributes(x$coefficients$a_domestic[[1]])$mcpar[3],
                                                           attributes(x$coefficients$a_domestic)$mcpar[3])})))
  
  mc_seq <- seq(from = mc_start, to = mc_end, by = mc_thin)[pos_thin]
  mc_start <- mc_seq[1]
  mc_thin <- mc_seq[2] - mc_seq[1]
  
  result$g0_i <- coda::mcmc(t(result$g0_i), start = mc_start, thin = mc_thin)
  result$g <- coda::mcmc(t(result$g), start = mc_start, thin = mc_thin)
  if (global) {
    result$h <- coda::mcmc(t(result$h), start = mc_start, thin = mc_thin)
  }
  if (deter) {
    result$d <- coda::mcmc(t(result$d), start = mc_start, thin = mc_thin)
  }
  result$sigma <- coda::mcmc(t(result$sigma), start = mc_start, thin = mc_thin)
  
  result$data$X <- NULL
  data_names <- NULL
  country_names <- NULL
  for (i in names(data)) {
    result$data$X <- cbind(result$data$X, data[[i]]$data$x_domestic)
    data_names <- c(data_names, dimnames(data[[i]]$data$x_domestic)[[2]])
    country_names <- c(country_names, rep(i, length(dimnames(data[[i]]$data$x_domestic)[[2]])))
  }
  dimnames(result$data$X)[[2]] <- data_names
  index <- data.frame("country" = country_names,
                      "variable" = data_names,
                      stringsAsFactors = FALSE)
  result$index <- index
  
  if (deter) {
    data_names <- NULL
    for (i in names(data)) {
      if (data[[i]]$model$type == "VEC" & !is.null(data[[i]]$data$deterministics)) {
        if (!is.null(data[[i]]$data$deterministics$unrestricted)) {
          result$data$deterministics <- cbind(result$data$deterministics, data[[i]]$data$deterministics$unrestricted)
          data_names <- c(data_names, dimnames(data[[i]]$data$deterministics$unrestricted)[[2]])
        }
        if (!is.null(data[[i]]$data$deterministics$restricted)) {
          result$data$deterministics <- cbind(result$data$deterministics, data[[i]]$data$deterministics$restricted) 
          data_names <- c(data_names, dimnames(data[[i]]$data$deterministics$restricted)[[2]])
        }
      } else {
        result$data$deterministics <- cbind(result$data$deterministics, data[[i]]$data$deterministics)
        data_names <- c(data_names, dimnames(data[[i]]$data$deterministics)[[2]])
      }
    }
    dimnames(result$data$deterministics)[[2]] <- data_names
    pos_det_name <- NULL
    pos_det <- NULL
    for (i in 1:length(data_names)) {
      if (!data_names[i] %in% pos_det_name) {
        pos_det_name <- c(pos_det_name, data_names[i])
        pos_det <- c(pos_det, i)
      }
    }
    temp_tsp <- stats::tsp(result$data$deterministics)
    result$data$deterministics <- stats::as.ts(as.matrix(result$data$deterministics[, pos_det]))
    dimnames(result$data$deterministics)[[2]] <- pos_det_name
    stats::tsp(result$data$deterministics) <- temp_tsp
  }
  
  if (global) {
    data_names <- NULL
    for (i in names(data)) {
      if (global_i[i]) {
        result$data$exogen <- cbind(result$data$exogen, data[[i]]$data$x_global)
        data_names <- c(data_names, dimnames(data[[i]]$data$x_global)[[2]])        
      }
    }
    dimnames(result$data$exogen)[[2]] <- data_names
    pos_det_name <- NULL
    pos_det <- NULL
    for (i in 1:length(data_names)) {
      if (!data_names[i] %in% pos_det_name) {
        pos_det_name <- c(pos_det_name, data_names[i])
        pos_det <- c(pos_det, i)
      }
    }
    temp_tsp <- stats::tsp(result$data$exogen)
    result$data$exogen <- stats::as.ts(as.matrix(result$data$exogen[, pos_det]))
    dimnames(result$data$exogen)[[2]] <- pos_det_name
    stats::tsp(result$data$exogen) <- temp_tsp
  }
  
  result$criteria <- criteria
  
  class(result) <- append("bgvar", class(result))
  return(result)
}