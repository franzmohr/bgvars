# Generate a VECX Model
# 
# Produces the model matrix of a VEC model with exogenous variables.
# 
# @param data containing the data and specification of a country model.
# 
# @return List containing the final model data and specifications.
# 
# @export
.gen_vecx <- function(object){
  
  k_domestic <- ncol(object$data$domestic)
  k_foreign <- ncol(object$data$foreign)
  r <- object$model$cointegration$rank
  
  level_domestic <- object$data$domestic
  level_foreign <- object$data$foreign
  ect <- cbind(level_domestic, level_foreign)
  names_ect <- c(paste(object$model$domestic$variables, ".l1", sep = ""),
                 paste("s.", object$model$foreign$variables, ".l1", sep = ""))
  
  global <- !is.null(object$data$global)
  if (global){
    p_global <- object$model$global$lags - 1
    k_global <- length(object$model$global$variables)
    level_global <- object$data$global
    ect <- cbind(ect, level_global)
    names_ect <- c(names_ect, paste(object$model$global$variables, ".l1", sep = "")) 
  }
  
  # Add restricted deterministic terms to error correction term
  det_res <- 0
  if (!is.null(object$data$deterministic$restricted)) {
    det_res <- ncol(object$data$deterministic$restricted)
    ect <- cbind(ect, object$data$deterministic$restricted)
    names_ect <- c(names_ect, dimnames(object$data$deterministic$restricted)[[2]])
  }
  k_ect <- dim(ect)[2]
  ect <- stats::lag(ect, -1)
  
  p_domestic <- object$model$domestic$lags - 1
  p_foreign <- object$model$foreign$lags - 1
  
  diff_domestic <- diff(level_domestic)
  total <- cbind(diff_domestic, ect)
  names_total <- c(paste("d.", object$model$domestic$variables, sep = ""), names_ect)
  
  # Lags of domestic variables
  n_d <- 0
  if (p_domestic > 0){
    for (i in 1:p_domestic){
      total <- cbind(total, stats::lag(diff_domestic, -i))
      names_total <- c(names_total, paste("d.", object$model$domestic$variables, ".l", i, sep = ""))
      n_d <- n_d + k_domestic
    }
  }
  
  # Lags of star variables
  diff_star <- diff(object$data$foreign)
  total <- cbind(total, diff_star)
  names_total <- c(names_total, paste("d.s.", object$model$foreign$variables, sep = ""))
  n_s <- k_foreign
  if (p_foreign > 0){
    for (i in 1:p_foreign){
      total <- cbind(total, stats::lag(diff_star, -i))
      names_total <- c(names_total, paste("d.s.", object$model$foreign$variables, ".l", i, sep = ""))
      n_s <- n_s + k_foreign
    }
  }
  
  # Lags of global variables
  n_g <- 0
  if (global){
    diff_global <- diff(object$data$global)
    total <- cbind(total, diff_global)
    names_total <- c(names_total, paste("d.", object$model$global$variables, sep = ""))
    n_g <- k_global
    if (p_global > 0){
      for (i in 1:p_global){
        total <- cbind(total, stats::lag(diff_global, -i))
        names_total <- c(names_total, paste("d.", object$model$global$variables, ".l", i, sep = ""))
        n_g <- n_g + k_global
      }
    }
  }
  
  # Add unrestrited deterministic terms
  det_unres <- 0
  if (!is.null(object$data$deterministic$unrestricted)) {
    det_unres <- ncol(object$data$deterministic$unrestricted)
    total <- cbind(total, object$data$deterministic$unrestricted)
    names_total <- c(names_total, dimnames(object$data$deterministic$unrestricted)[[2]])
  }
  
  total <- stats::na.omit(total)
  dimnames(total)[[2]] <- names_total
  tsp_series <- stats::tsp(total)
  time_series <- stats::time(total)
  
  # Final data preparations
  used_t <- seq(to = nrow(total), length.out = nrow(object$data$domestic) - object$model$max_lag)
  ts_start = time_series[used_t][1]
  y <- stats::ts(as.matrix(total[used_t, 1:k_domestic]),
                 start = ts_start, frequency = tsp_series[3],
                 class = c("mts", "ts", "matrix"))
  dimnames(y) <- list(NULL, names_total[1:k_domestic])
  ect <- stats::ts(as.matrix(total[used_t, k_domestic + 1:k_ect]),
                 start = ts_start, frequency = tsp_series[3],
                 class = c("mts", "ts", "matrix"))
  x <- stats::ts(as.matrix(total[used_t, -(1:(k_domestic + k_ect))]),
                 start = ts_start, frequency = tsp_series[3],
                 class = c("mts", "ts", "matrix"))
  
  sur <- kronecker(cbind(ect, x), diag(1, k_domestic))
  
  structural <- object$model[["structural"]]
  if (structural & k_domestic > 1) {
    y_structural <- kronecker(-y, diag(1, k_domestic))
    pos <- NULL
    for (j in 1:k_domestic) {
      pos <- c(pos, (j - 1) * k_domestic + 1:j)
    }
    sur <- cbind(sur, y_structural[, -pos])
  }
  
  result <- list(Y = y, W = ect, X = x, SUR = sur)
  return(result)
}