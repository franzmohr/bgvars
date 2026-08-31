
.create_weights_submodel <- function(submodel, period = NULL) {
  
  endogen <- submodel[["endogen"]]
  weights <- submodel[["weights"]]
  
  use_rolling_window <- length(period) == 1
  
  if ("ts" %in% class(weights)) {
    
    # Time stamps of endogenous data
    tt <- NROW(endogen)
    endogen_time <- as.numeric(stats::time(endogen))
    endogen_tsp <- stats::tsp(endogen)
    
    # Time stamps of weight data
    weight_time <- as.numeric(stats::time(weights))
    
    # Determine availability of weight data per period
    if (use_rolling_window) {
      availability_matrix <- matrix(NA, tt, period)
      for (i in 1:tt) {
        if (endogen_time[i] <= weight_time[period]) {
          availability_matrix[i,] <- weight_time[1:period]
        }
        if (endogen_time[i] >= weight_time[period]) {
          if (any(floor(endogen_time[i]) == weight_time)) {
            pos_t <- which(floor(endogen_time[i]) == weight_time)
            pos_t <- (pos_t - period + 1):pos_t
          }
          # If condition is not met, the pos_t, from the last iteration will be used.
          availability_matrix[i,] <- weight_time[pos_t]
        } 
      } 
    } else {
      availability_matrix <- t(matrix(period, nrow = length(period), ncol = tt))
    }
    
    # Create a weight matrix for each period
    w <- stats::ts(matrix(NA, tt, ncol(weights)),
                   start = endogen_tsp[1], frequency = endogen_tsp[3])
    dimnames(w)[[2]] <- dimnames(weights)[[2]]
    
    for (i in 1:tt){
      temp <- colSums(weights[stats::time(weights) %in% as.character(availability_matrix[i,]), ])
      w[i, ] <- temp / sum(temp)
    }
    
  } else {
    stop("Element 'weights' must be of class 'ts'.")
  }
  
  return(w)
  
}
