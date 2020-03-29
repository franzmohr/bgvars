# Generate a VECX Model
# 
# Produces the model matrix of a VEC model with exogenous variables.
# 
# @param data containing the data and specification of a country model.
# 
# @return List containing the final model data and specifications.
# 
# @export
.gen_vecx <- function(data){
  k_domestic <- ncol(data$data$x_domestic)
  k_foreign <- ncol(data$data$x_foreign)
  r <- data$model$cointegration$rank
  
  level_domestic <- data$data$x_domestic
  level_foreign <- data$data$x_foreign
  ect <- cbind(level_domestic, level_foreign)
  names_ect <- c(paste(data$model$domestic$variables, ".l1", sep = ""),
                 paste("s.", data$model$foreign$variables, ".l1", sep = ""))
  
  global <- !is.null(data$data$x_global)
  if (global){
    p_global <- data$model$global$lags - 1
    k_global <- length(data$model$global$variables)
    level_global <- data$data$x_global
    ect <- cbind(ect, level_global)
    names_ect <- c(names_ect, paste(data$model$global$variables, ".l1", sep = "")) 
  }
  
  # Add restricted deterministic terms to error correction term
  det_res <- 0
  if (!is.null(data$data$deterministics$restricted)) {
    det_res <- ncol(data$data$deterministics$restricted)
    ect <- cbind(ect, data$data$deterministics$restricted)
    names_ect <- c(names_ect, dimnames(data$data$deterministics$restricted)[[2]])
  }
  k_ect <- dim(ect)[2]
  ect <- stats::lag(ect, -1)
  
  p_domestic <- data$model$domestic$lags - 1
  p_foreign <- data$model$foreign$lags - 1
  
  diff_domestic <- diff(level_domestic)
  total <- cbind(diff_domestic, ect)
  names_total <- c(paste("d.", data$model$domestic$variables, sep = ""), names_ect)
  
  # Lags of domestic variables
  n_d <- 0
  if (p_domestic > 0){
    for (i in 1:p_domestic){
      total <- cbind(total, stats::lag(diff_domestic, -i))
      names_total <- c(names_total, paste("d.", data$model$domestic$variables, ".l", i, sep = ""))
      n_d <- n_d + k_domestic
    }
  }
  
  # Lags of star variables
  diff_star <- diff(data$data$x_foreign)
  total <- cbind(total, diff_star)
  names_total <- c(names_total, paste("d.s.", data$model$foreign$variables, sep = ""))
  n_s <- k_foreign
  if (p_foreign > 0){
    for (i in 1:p_foreign){
      total <- cbind(total, stats::lag(diff_star, -i))
      names_total <- c(names_total, paste("d.s.", data$model$foreign$variables, ".l", i, sep = ""))
      n_s <- n_s + k_foreign
    }
  }
  
  # Lags of global variables
  n_g <- 0
  if (global){
    diff_global <- diff(data$data$x_global)
    total <- cbind(total, diff_global)
    names_total <- c(names_total, paste("d.", data$model$global$variables, sep = ""))
    n_g <- k_global
    if (p_global > 0){
      for (i in 1:p_global){
        total <- cbind(total, stats::lag(diff_global, -i))
        names_total <- c(names_total, paste("d.", data$model$global$variables, ".l", i, sep = ""))
        n_g <- n_g + k_global
      }
    }
  }
  
  # Add unrestrited deterministic terms
  det_unres <- 0
  if (!is.null(data$data$deterministics$unrestricted)) {
    det_unres <- ncol(data$data$deterministics$unrestricted)
    total <- cbind(total, data$data$deterministics$unrestricted)
    names_total <- c(names_total, dimnames(data$data$deterministics$unrestricted)[[2]])
  }
  
  total <- stats::na.omit(total)
  dimnames(total)[[2]] <- names_total
  
  # Final data preparations
  used_t <- seq(to = nrow(total), length.out = nrow(data$data$x_domestic) - data$model$max_lag)
  y <- t(total[used_t, 1:k_domestic])
  dimnames(y) <- list(names_total[1:k_domestic], NULL)
  ect <- t(total[used_t, k_domestic + 1:(k_ect)])
  x <- t(total[used_t, -(1:(k_domestic + k_ect))])
  
  result <- list(y = y, ect = ect, x = x)
  result$domestic <- list("dim" = k_domestic, "lag" = p_domestic)
  result$foreign <- list("dim" = k_foreign, "lag" = p_foreign)
  if (global) {
    result$global <- list("dim" = k_global, "lag" = p_global)
  }
  result$deterministics <- list("restricted" = det_res, "unresticted" = det_unres)
  return(result)
}