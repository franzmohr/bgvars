# Generate a VARX Model
# 
# Produces the model matrix of a VEC model with exogenous variables.
# 
# @param data containing the data and specification of a country model.
# 
# @return List containing the final model data and specifications.
# 
# @export
.gen_varx <- function(object){
  
  # rm(list = ls()[-which(ls() == "object")])
  
  k_domestic <- ncol(object$data$domestic)
  k_foreign <- ncol(object$data$foreign)
  
  level_domestic <- object$data$domestic
  level_foreign <- object$data$foreign
  
  global <- !is.null(object$data$global)
  if (global){
    p_global <- object$model$global$lags
    k_global <- length(object$model$global$variables)
    level_global <- object$data$global
  }
  
  p_domestic <- object$model$domestic$lags
  p_foreign <- object$model$foreign$lags
  
  total <- level_domestic
  names_total <- object$model$domestic$variables
  
  # Lags of domestic variables
  n_d <- 0
  if (p_domestic > 0){
    for (i in 1:p_domestic){
      total <- cbind(total, stats::lag(level_domestic, -i))
      names_total <- c(names_total, paste(object$model$domestic$variables, ".l", i, sep = ""))
      n_d <- n_d + k_domestic
    }
  }
  
  # Lags of star variables
  total <- cbind(total, level_foreign)
  names_total <- c(names_total, paste("s.", object$model$foreign$variables, sep = ""))
  n_s <- k_foreign
  if (p_foreign > 0){
    for (i in 1:p_foreign){
      total <- cbind(total, stats::lag(level_foreign, -i))
      names_total <- c(names_total, paste("s.", object$model$foreign$variables, ".l", i, sep = ""))
      n_s <- n_s + k_foreign
    }
  }
  
  # Lags of global variables
  n_g <- 0
  if (global){
    total <- cbind(total, level_global)
    names_total <- c(names_total, paste(object$model$global$variables, sep = ""))
    n_g <- k_global
    if (p_global > 0){
      for (i in 1:p_global){
        total <- cbind(total, stats::lag(level_global, -i))
        names_total <- c(names_total, paste(object$model$global$variables, ".l", i, sep = ""))
        n_g <- n_g + k_global
      }
    }
  }
  
  # Add deterministic terms
  k_det <- 0
  if (!is.null(object$data$deterministic)) {
    temp <- object$data$deterministic
    name_temp <- dimnames(object$data$deterministic)[[2]]
    k_det <- ncol(temp)
    total <- cbind(total, temp)
    names_total <- c(names_total, name_temp)
    rm(list = c("temp", "name_temp"))
  }
  
  total <- stats::na.omit(total)
  dimnames(total)[[2]] <- names_total
  tsp_series <- stats::tsp(total)
  time_series <- stats::time(total)
  
  # Final object preparations
  used_t <- seq(to = nrow(total), length.out = nrow(object$data$domestic) - object$model$max_lag)
  ts_start = time_series[used_t][1]
  y <- stats::ts(as.matrix(total[used_t, 1:k_domestic]),
                 start = ts_start, frequency = tsp_series[3],
                 class = c("mts", "ts", "matrix"))
  dimnames(y) <- list(NULL, names_total[1:k_domestic])
  x <- stats::ts(as.matrix(total[used_t, -(1:k_domestic)]),
                 start = ts_start, frequency = tsp_series[3],
                 class = c("mts", "ts", "matrix"))
  sur <- kronecker(x, diag(1, k_domestic))
  
  # Add structural variables to SUR representation here
  structural <- object$model$structural
  if (structural & k_domestic > 1) {
    y_structural <- kronecker(-y, diag(1, k_domestic))
    pos <- NULL
    for (j in 1:k_domestic) {
      pos <- c(pos, (j - 1) * k_domestic + 1:j)
    }
    sur <- cbind(sur, y_structural[, -pos])
  }

  return(list(Y = y, Z = x, SUR = sur))
}