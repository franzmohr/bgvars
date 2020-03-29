#' Differences of Variables
#' 
#' Calculates the first difference of variables in a time-series object or in a
#' list of named time-series objects
#' 
#' @param data a named time-series object or a list of named time-series objects.
#' @param variables a character vector of variables that should be differenced.
#' If \code{NULL} (default), all variables are differenced.
#' 
#' @return A differenced time-series object or a named list of differenced time-series objects.
#' 
#' @examples 
#' data("gvar2016")
#' 
#' country_data <- gvar2016$country_data
#' 
#' # Take first difference of all variables y and Dp
#' country_data <- diff_variables(country_data, variables = c("y", "Dp"))
#' 
#' @export
diff_variables <- function(data, variables = NULL){
  if (any(class(data) == "list")) {
    if (is.null(names(data))) {
      stop("Argument 'data' must be named, if a list is provided.")
    }
    if(sum(unlist(lapply(data, class)) == "ts") != length(data)) {stop("Data must be of class 'ts'.")} 
    data <- lapply(data, .diff_func, variables)
  } else {
    if(!"ts" %in% class(data)) {
      stop("Data must be of class 'ts'.")
    } else {
      data <- .diff_func(data, variables) 
    }
  }
  
  return(data)
}

.diff_func <- function(x, variables){
  tsp_all <- stats::tsp(x)
  tsp_all[1] <- tsp_all[1] + 1 / tsp_all[3]
  if (is.null(variables)){
    for (i in dimnames(x)[[2]]){
      x[, i] <- c(NA, diff(x[, i]))
    }
  } else {
    for (i in variables){
      if (is.element(i, dimnames(x)[[2]])){
        x[, i] <- c(NA, diff(x[, i]))
      }
    } 
  }
  x <- x[-1, ]
  x <- stats::ts(x, start = tsp_all[1], frequency = tsp_all[3])
  return(x)
}