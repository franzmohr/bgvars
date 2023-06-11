#' Create Regional Series
#' 
#' Combines multiple country-specific time series to regional series.
#' 
#' @param country_data a named list of time-series objects containing country data.
#' @param weight_data a named list of named time-series objects. See 'Details'.
#' @param region_weights a multivariate time-series object containing data used
#' to weight the observations in \code{country_data}.
#' @param regions a named list of character vectors containing specifications for countries
#' in \code{country_data}, which should be combined to a region. The name of the
#' respective list element will become the name of the region.
#' @param period either a single integer or a numeric vector specifiying the periods in
#' \code{region_weights}, which should be used to construct weights.
#' 
#' @details
#' 
#' If a numeric vector is provided for argument \code{period}, the function weights
#' country-specific observations based on the sums over the specified periods.
#' If an integer is proved, the country-specific observations are weighted
#' according to rolling sums over the last \code{period} periods. If country
#' data start earlier than the series in \code{region_weights}, the sums over
#' the first \code{period} observations of \code{region_weights} are used
#' until the periods match.
#' 
#' The elements of \code{weight_data} should be time-series objects,
#' where colums are named after the trading partners of the respective country
#' and the names of the elements should match the country names in
#' \code{country_data}. For each country list element the time-series object should
#' contain zeros in the column of the respective country. For examle, the column
#' "US" in the "US"-element of the list should only contain zeros.
#' 
#' @return A list containing the following elements:
#' \item{country_data}{a named list of updated time-series objects with country data.}
#' \item{weight_data}{a named list of named time-series objects for the
#' generation of country-specific weight matrices.}
#' 
#' @examples
#' # Load data
#' data("gvar2016")
#' 
#' # Create regions
#' temp <- create_regions(country_data = gvar2016$country_data,
#'              weight_data = gvar2016$weight_data,
#'              region_weights = gvar2016$region_weights,
#'              regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")),
#'              period = 3)
#' 
#' # New data sets
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data 
#' 
#' @export
create_regions <- function(country_data, weight_data, region_weights, regions, period){
  
  tt <- unique(unlist(lapply(country_data, NROW)))
  if (length(tt) > 1) {stop("Country data must have the same numbers of observations.")}
  
  if ((!"list" %in% class(regions)) | is.null(names(regions))) {stop("Object 'regions' must be a named list.")}
  
  if (length(unique(unlist(regions))) < length(unlist(regions))) {
    stop("The same country is not allowed to be in more than one region.") 
  }
  
  vars <- unique(unlist(lapply(country_data, function(x){return(dimnames(x)[[2]])})))
  
  var_exist <- matrix(FALSE, length(country_data), length(vars))
  dimnames(var_exist) <- list(names(country_data), vars)
  for (i in names(country_data)) {
    var_exist[i, dimnames(country_data[[i]])[[2]]] <- TRUE
  }
  
  if (length(period) == 1) {
    t.temp <- as.numeric(stats::time(country_data[[1]]))
    t.avail <- as.numeric(stats::time(region_weights))
    
    t <- matrix(NA, tt, period)
    for (i in 1:tt) {
      if (t.temp[i] <= t.avail[period]) {
        t[i,] <- t.avail[1:period]
      }
      # Use last available values of region weights if country data is more recent
      if (t.temp[i] >= t.avail[period]) {
        if (any(floor(t.temp[i]) == t.avail)) {
          pos_t <- which(floor(t.temp[i]) == t.avail)
          pos_t <- (pos_t - period + 1):pos_t 
        }
        # If condition is not met, the pos_t, from the last iteration will be used.
        t[i,] <- t.avail[pos_t]
      }
    }
  }
  
  r_names <- names(regions)
  all_r_countries <- unlist(regions)
  names(all_r_countries) <- NULL
  r_tsp <- stats::tsp(country_data[[1]])
  
  r.temp <- c()
  for (i in 1:length(regions)) {
    vars_r <- apply(var_exist[regions[[i]], ], 2, any)
    vars_r <- dimnames(var_exist)[[2]][vars_r]
    
    r_temp <- stats::ts(matrix(NA, tt, length(vars_r)), start = r_tsp[1], frequency = r_tsp[3])
    dimnames(r_temp)[[2]] <- vars_r
    
    for (j in vars_r) {
      c_temp <- matrix(NA, tt, length(regions[[i]]))
      dimnames(c_temp)[[2]] <- regions[[i]]
      for (k in regions[[i]]) {
        if (var_exist[k , j]) {
          c_temp[, k] <- country_data[[k]][, j] 
        }
      }
      c_temp <- c_temp[, var_exist[regions[[i]], j]]
      
      if (NCOL(c_temp) > 1) {
        if (length(period) == 1) {
          for (k in 1:tt) {
            temp <- colSums(region_weights[which(dimnames(region_weights)[[1]] %in% t[k,]), dimnames(c_temp)[[2]]])
            temp <- temp / sum(temp)
            r_temp[k, j] <- sum(c_temp[k, ] * temp)
          }
        } else {
          temp <- colSums(region_weights[which(dimnames(region_weights)[[1]] %in% as.character(period)), dimnames(c_temp)[[2]]])
          temp <- temp / sum(temp)
          r_temp[, j] <- c_temp %*% matrix(temp)
        }
      } else {
        r_temp[, j] <- c_temp
      }
    }
    r.temp <- c(r.temp, list(r_temp))
  }
  names(r.temp) <- names(regions)
  
  w_temp <- c()
  for (i in names(country_data)) {
    if (!i %in% all_r_countries) {
      w_i <- weight_data[[i]][, -which(dimnames(weight_data[[i]])[[2]] %in% all_r_countries)]
      w_i_names <- c(dimnames(w_i)[[2]], names(regions))
      w_i <- cbind(w_i, matrix(0, ncol = length(regions)))
      dimnames(w_i)[[2]] <- w_i_names
      
      for (j in r_names) {
        w_i[, j] <- rowSums(weight_data[[i]][, regions[[j]]])
      }
      w_temp <- c(w_temp, list(w_i))
      rm(w_i)
    }
  }
  
  tot_names <- unique(unlist(lapply(weight_data, function(x){dimnames(x)[[2]]})))
  tot_names <- tot_names[-which(tot_names %in% all_r_countries)]
  for (i in r_names) {
    w_i <- weight_data[[regions[[i]][1]]][, tot_names] * 0
    w_i <- cbind(w_i, matrix(0, nrow(w_i), length(regions)))
    dimnames(w_i)[[2]] <- c(tot_names, r_names)
    for (j in regions[[i]]) {
      # Add non-regional data
      w_i[, tot_names] <- w_i[, tot_names] + weight_data[[j]][, tot_names]
      # Add regional data
      for (k in r_names) {
        if (k != i) {
          for (l in regions[[k]]) {
            w_i[, k] <- w_i[, k] + weight_data[[j]][, l] 
          }
        }
      }
    }
    w_temp <- c(w_temp, list(w_i)) 
  }
  
  data <- c()
  data.names <- c()
  for (i in names(country_data)) {
    if (!i %in% all_r_countries) {
      data <- c(data , list(country_data[[i]]))
      data.names <- c(data.names, i)
    }
  }
  for (i in names(r.temp)) {
    data <- c(data, list(r.temp[[i]]))
    data.names <- c(data.names, i)
  }
  names(data) <- data.names
  names(w_temp) <- data.names
  
  result <- list("country_data" = data,
                 "weight_data" = w_temp)
  return(result)
}