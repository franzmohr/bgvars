#' Create Regional Series
#' 
#' Combines multiple country-specific time series to regional series.
#' 
#' @param submodel_data a named list of class 'submodeldata'.
#' @param region_weights a multivariate time-series object containing data used
#' to weight the observations in \code{submodel_data}.
#' @param regions a named list of character vectors containing specifications for countries
#' in \code{submodel_data}, which should be combined to a region. The name of the
#' respective list element will become the name of the region.
#' @param period either a single integer or a numeric vector specifiying the periods in
#' \code{region_weights}, which should be used to construct weights.
#' 
#' @details
#' 
#' If a numeric vector is provided for argument \code{period}, the function weights
#' submodel-specific observations based on the sums over the specified periods.
#' If an integer is proved, the submodel-specific observations are weighted
#' according to rolling sums over the last \code{period} periods. If submodel
#' data starts earlier than the series in \code{region_weights}, the sums over
#' the first \code{period} observations of \code{region_weights} are used
#' until the periods match.
#' 
#' @return An object of class 'submodeldata'.
#' 
#' @examples
#' # Load data
#' data("gvar2019")
#' 
#' new_regions <- list("EA" = c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL"))
#' 
#' # Create regions
#' submodel_data <- create_regions(submodel_data = gvar2019[["submodel_data"]],
#'                                 region_weights = gvar2019[["region_weights"]],
#'                                 regions = new_regions,
#'                                 period = 3)
#' 
#' @export
create_regions <- function(submodel_data, region_weights, regions, period){
  
  tt <- unique(unlist(lapply(submodel_data, function(x) {NROW(x[["endogen"]])})))
  if (length(tt) > 1) {stop("Currently, submodel data must have the same numbers of observations.")}
  
  if ((!"list" %in% class(regions)) | is.null(names(regions))) {stop("Object 'regions' must be a named list.")}
  
  if (length(unique(unlist(regions))) < length(unlist(regions))) {
    stop("The same country is not allowed to be in more than one region.") 
  }
  
  vars <- unique(unlist(lapply(submodel_data, function(x){return(dimnames(x[["endogen"]])[[2]])})))
  
  # Create a table containing information on which variable is availabe within
  # a submodel.
  submodel_names <- names(submodel_data)
  var_exist <- matrix(FALSE, length(submodel_data), length(vars))
  dimnames(var_exist) <- list(submodel_names, vars)
  for (i in submodel_names) {
    var_exist[i, dimnames(submodel_data[[i]][["endogen"]])[[2]]] <- TRUE
  }
  
  # Create a matrix, where each row contains the oberservations that should be
  # used for weight construction for each period
  if (length(period) == 1) {
    t_temp <- as.numeric(stats::time(submodel_data[[1]][["endogen"]]))
    t_avail <- as.numeric(stats::time(region_weights))
    
    rolling_window_weights <- matrix(NA, tt, period)
    for (i in 1:tt) {
      if (t_temp[i] <= t_avail[period]) {
        rolling_window_weights[i,] <- t_avail[1:period]
      }
      # Use last available values of region weights if country data is more recent
      if (t_temp[i] >= t_avail[period]) {
        if (any(floor(t_temp[i]) == t_avail)) {
          pos_t <- which(floor(t_temp[i]) == t_avail)
          pos_t <- (pos_t - period + 1):pos_t 
        }
        # If condition is not met, the pos_t, from the last iteration will be used.
        rolling_window_weights[i,] <- t_avail[pos_t]
      }
    }
  }
  
  r_names <- names(regions)
  all_r_countries <- unlist(regions)
  names(all_r_countries) <- NULL
  r_tsp <- stats::tsp(submodel_data[[1]][["endogen"]])
  
  # Create submodel entries for regions, but do not add them to main data set yet.
  endogen_temp <- c()
  for (i in 1:length(regions)) {
    
    # Check variable availability
    vars_r <- apply(var_exist[regions[[i]], ], 2, any)
    vars_r <- dimnames(var_exist)[[2]][vars_r]
    
    r_temp <- stats::ts(matrix(NA, tt, length(vars_r)), start = r_tsp[1], frequency = r_tsp[3])
    dimnames(r_temp)[[2]] <- vars_r
    
    for (j in vars_r) {
      # Create matrix for an individual series with the data from all submodels,
      # from which data is used to construct the group
      c_temp <- matrix(NA, tt, length(regions[[i]]))
      dimnames(c_temp)[[2]] <- regions[[i]]
      for (k in regions[[i]]) {
        if (var_exist[k , j]) {
          c_temp[, k] <- submodel_data[[k]][["endogen"]][, j] 
        }
      }
      # Only use submodels for which observations are available for that variable
      c_temp <- c_temp[, var_exist[regions[[i]], j]]
      
      if (NCOL(c_temp) > 1) {
        # In case of rolling window weights calculate each observation separately
        if (length(period) == 1) {
          for (k in 1:tt) {
            # Create weights
            temp <- colSums(region_weights[which(dimnames(region_weights)[[1]] %in% rolling_window_weights[k,]), dimnames(c_temp)[[2]]])
            temp <- temp / sum(temp)
            # Calculate weighted mean
            r_temp[k, j] <- sum(c_temp[k, ] * temp)
          }
        } else {
          # Create weights
          temp <- colSums(region_weights[which(dimnames(region_weights)[[1]] %in% as.character(period)), dimnames(c_temp)[[2]]])
          temp <- temp / sum(temp)
          r_temp[, j] <- c_temp %*% matrix(temp)
        }
      } else {
        r_temp[, j] <- c_temp
      }
    }
    endogen_temp <- c(endogen_temp, list(r_temp))
  }
  names(endogen_temp) <- names(regions)
  
  # Update weight data
  w_temp <- c()
  for (i in names(submodel_data)) {
    if (!i %in% all_r_countries) {
      w_i <- submodel_data[[i]][["weights"]][, -which(dimnames(submodel_data[[i]][["weights"]])[[2]] %in% all_r_countries)]
      w_i_names <- c(dimnames(w_i)[[2]], names(regions))
      w_i <- cbind(w_i, matrix(0, ncol = length(regions)))
      dimnames(w_i)[[2]] <- w_i_names
      
      for (j in r_names) {
        w_i[, j] <- rowSums(submodel_data[[i]][["weights"]][, regions[[j]]])
      }
      w_temp <- c(w_temp, list(w_i))
      rm(w_i)
    }
  }
  
  tot_names <- unique(unlist(lapply(submodel_data, function(x){dimnames(x[["weights"]])[[2]]})))
  tot_names <- tot_names[-which(tot_names %in% all_r_countries)]
  for (i in r_names) {
    w_i <- submodel_data[[regions[[i]][1]]][["weights"]][, tot_names] * 0
    w_i <- cbind(w_i, matrix(0, nrow(w_i), length(regions)))
    dimnames(w_i)[[2]] <- c(tot_names, r_names)
    for (j in regions[[i]]) {
      # Add non-regional data
      w_i[, tot_names] <- w_i[, tot_names] + submodel_data[[j]][["weights"]][, tot_names]
      # Add regional data
      for (k in r_names) {
        if (k != i) {
          for (l in regions[[k]]) {
            w_i[, k] <- w_i[, k] + submodel_data[[j]][["weights"]][, l] 
          }
        }
      }
    }
    w_temp <- c(w_temp, list(w_i)) 
  }
  
  data <- c()
  data.names <- c()
  for (i in names(submodel_data)) {
    if (!i %in% all_r_countries) {
      data <- c(data , list(submodel_data[[i]][["endogen"]]))
      data.names <- c(data.names, i)
    }
  }
  for (i in names(endogen_temp)) {
    data <- c(data, list(endogen_temp[[i]]))
    data.names <- c(data.names, i)
  }
  names(data) <- data.names
  names(w_temp) <- data.names
  
  result <- NULL
  for (i in data.names) {
    result[[i]] <- list("endogen" = data[[i]],
                        "weights" = w_temp[[i]])
  }
  class(result) <- list("submodeldata", "list")
  
  return(result)
}