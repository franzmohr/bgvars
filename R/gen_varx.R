# Generate a VARX Model
# 
# Produces the model matrix of a VEC model with exogenous variables.
# 
# @param data containing the data and specification of a country model.
# 
# @return List containing the final model data and specifications.
# 
# @export
.gen_varx <- function(data){
  k_domestic <- ncol(data$data$x_domestic)
  k_foreign <- ncol(data$data$x_foreign)
  
  level_domestic <- data$data$x_domestic
  level_foreign <- data$data$x_foreign
  
  global <- !is.null(data$data$x_global)
  if (global){
    p_global <- data$model$global$lags
    k_global <- length(data$model$global$variables)
    level_global <- data$data$x_global
  }
  
  p_domestic <- data$model$domestic$lags
  p_foreign <- data$model$foreign$lags
  
  total <- level_domestic
  names_total <- data$model$domestic$variables
  
  # Lags of domestic variables
  n_d <- 0
  if (p_domestic > 0){
    for (i in 1:p_domestic){
      total <- cbind(total, stats::lag(level_domestic, -i))
      names_total <- c(names_total, paste(data$model$domestic$variables, ".l", i, sep = ""))
      n_d <- n_d + k_domestic
    }
  }
  
  # Lags of star variables
  total <- cbind(total, level_foreign)
  names_total <- c(names_total, paste("s.", data$model$foreign$variables, sep = ""))
  n_s <- k_foreign
  if (p_foreign > 0){
    for (i in 1:p_foreign){
      total <- cbind(total, stats::lag(level_foreign, -i))
      names_total <- c(names_total, paste("s.", data$model$foreign$variables, ".l", i, sep = ""))
      n_s <- n_s + k_foreign
    }
  }
  
  # Lags of global variables
  n_g <- 0
  if (global){
    total <- cbind(total, level_global)
    names_total <- c(names_total, paste(data$model$global$variables, sep = ""))
    n_g <- k_global
    if (p_global > 0){
      for (i in 1:p_global){
        total <- cbind(total, stats::lag(level_global, -i))
        names_total <- c(names_total, paste(data$model$global$variables, ".l", i, sep = ""))
        n_g <- n_g + k_global
      }
    }
  }
  
  # Add unrestrited deterministic terms
  k_det <- 0
  if (!is.null(data$data$deterministics)) {
    temp <- data$data$deterministics
    name_temp <- dimnames(data$data$deterministics)[[2]]
    k_det <- ncol(temp)
    total <- cbind(total, temp)
    names_total <- c(names_total, name_temp)
    rm(list = c("temp", "name_temp"))
  }
  
  total <- stats::na.omit(total)
  dimnames(total)[[2]] <- names_total
  
  # Final data preparations
  used_t <- seq(to = nrow(total), length.out = nrow(data$data$x_domestic) - data$model$max_lag)
  y <- t(total[used_t, 1:k_domestic])
  dimnames(y) <- list(names_total[1:k_domestic], NULL)
  x <- t(total[used_t, -(1:k_domestic)])
  
  result <- list(y = y, x = x)
  result$domestic <- list("dim" = k_domestic, "lag" = p_domestic)
  result$foreign <- list("dim" = k_foreign, "lag" = p_foreign)
  if (global) {
    result$global <- list("dim" = k_global, "lag" = p_global)
  }
  
  result$deterministics <- list("dim" = k_det) 
  return(result)
}