#' Create Weight Matrices
#' 
#' Produces a matrix of fixed or an array of time varying country-specific
#' weights for a GVAR model.
#' 
#' @param weight_data a named list of named time-series objects or a named
#' matrix containing data for the construction of weights. See 'Details'.
#' @param period either an integer for time varying weights or a numeric
#' vector specifying the number of past periods in \code{weight_data} that
#' should be used to calculate constant weights. See 'Details'.
#' @param country_data a named list of time-series objects containing country
#' data. Only requried if \code{period} is an integer.
#' 
#' @details 
#' 
#' The function assists in the creation of country-specific weight matrices.
#' If a numeric vector is provided as \code{period}, the function calculates
#' country weights based on the sums over the specified periods. If an integer
#' is proved, the weights are constructed from rolling sums over the last
#' \code{period} periods. If a country series begins earlier than its
#' corresponding weight series, the sums over the first \code{period}
#' observations of \code{weight_data} are used until the periods match.
#' 
#' If \code{weight_data} is a list, its elements should be time-series objects,
#' where colums should be named after the trading partner. The names of the list
#' elements should match the country names in \code{country_data}. If a matrix is
#' provided, its rows should contain the name of the reference country and the
#' columns should be named after the trading partner. All cells should contain
#' zeros when the reference country and the trading partner are the same,
#' i.e. row and column names are equal.
#' 
#' @return A named matrix or an array containing country-specific weight matrices.
#' 
#' @examples
#' # Load data
#' data("gvar2019")
#'
#' # Data objects
#' country_data <- gvar2019[["country_data"]]
#' weight_data <- gvar2019[["weight_data"]]
#' 
#' weight_data <- create_weights(weight_data = weight_data, period = 3,
#'                               country_data = country_data)
#' 
#' @export
create_weights <- function(weight_data, period, country_data = NULL){
  # Check if weight.data is a list or matrix.
  if (any(class(weight_data) %in% c("matrix", "list"))) {
    if ("list" %in% class(weight_data)) {
      l <- TRUE
    } else {
      l <- FALSE
    } 
  } else {
    stop("Object 'weight_data' must be of class list or matrix.")
  }
  
  if (l) {
    if (length(period) > 1) {
      period <- as.character(period)
      if (l){
        w <- matrix(NA, length(weight_data), length(weight_data))
        dimnames(w) <- list(names(weight_data), names(weight_data))
        for (i in names(weight_data)) {
          temp <- colSums(weight_data[[i]][stats::time(weight_data[[i]]) %in% period, ])
          w[i, names(temp)] <- temp / sum(temp)
        }
      }
    } else {
      if (is.null(country_data)){
        stop("The construction of time varying weights requries to specify the 'country_data' argument.")
      }
      tt <- unique(unlist(lapply(country_data, NROW)))
      if (length(tt) > 1) {
        stop("Objects in 'country_data' do not have the same number of observations.")
      }
      t_temp <- as.numeric(stats::time(country_data[[1]]))
      
      t_ind <- stats::tsp(country_data[[1]])
      t_avail <- unique(unlist(lapply(weight_data, NROW)))
      if (length(t_avail) > 1) {
        stop("Objects in 'weight_data' do not have the same number of observations.")
      }
      t_avail <- as.numeric(stats::time(weight_data[[1]]))
      
      t <- matrix(NA, tt, period)
      for (i in 1:tt) {
        if (t_temp[i] <= t_avail[period]) {
          t[i,] <- t_avail[1:period]
        }
        if (t_temp[i] >= t_avail[period]) {
          if (any(floor(t_temp[i]) == t_avail)) {
            pos_t <- which(floor(t_temp[i]) == t_avail)
            pos_t <- (pos_t - period + 1):pos_t
          }
          # If condition is not met, the pos_t, from the last iteration will be used.
          t[i,] <- t_avail[pos_t]
        }
      }
      
      w <- array(NA, dim = c(length(weight_data), length(weight_data), tt))
      dimnames(w) <- list(names(weight_data), names(weight_data), as.character(t_temp))
      for (i in 1:tt){
        for (j in names(country_data)) {
          temp <- colSums(weight_data[[j]][stats::time(weight_data[[j]]) %in% as.character(t[i,]), ])
          w[j, names(temp), i] <- temp / sum(temp)
        }
      }
    } 
  } else {
    if (is.null(dimnames(weight_data)[[1]]) | is.null(dimnames(weight_data)[[2]])) {
      stop("Both the rows and columns of the 'weight_data' must be named.")
    }
    
    r_names <- dimnames(weight_data)[[1]]
    r_names <- r_names[order(r_names)]
    c_names <- dimnames(weight_data)[[2]]
    c_names <- c_names[order(c_names)]
    weight_data <- weight_data[r_names, c_names]
    
    if (any(diag(weight_data) != 0)) {
      stop("Entries, where row and column names are the same must contain the value zero.")
    }
    
    w <- matrix(NA, NROW(weight_data), NCOL(weight_data))
    dimnames(w) <- list(dimnames(weight_data)[[1]], dimnames(weight_data)[[1]])
    for (i in dimnames(weight_data)[[1]]) {
      temp <- weight_data[i, ]
      w[i, names(temp)] <- temp / sum(temp)
    }
  }
  return(w)
}
