#' Differences of Variables
#' 
#' Produces the first difference of variables in an object of class 'submodeldata'.
#' 
#' @param data an object of class 'submodeldata'.
#' @param variables a character vector of variables that should be differenced, if
#' they appear in a time-series object. If \code{NULL} (default), all variables are differenced.
#' @param multi optional. Numeric by which the differenced series should be multiplicated.
#' 
#' @return A differenced time-series object or a list of differenced time-series objects.
#' 
#' @examples 
#' # Load data
#' data("gvar2019")
#' submodel_data <- gvar2019$submodel_data
#' 
#' # Take first difference of the variables "y" and "Dp" across all
#' # elements of object "submodel_data" and multiply them by 100
#' submodel_data <- diff(submodel_data, variables = c("y", "Dp"), multi = 100)
#' 
#' @export
diff.submodeldata <- function(data, variables = NULL, multi = NULL){
  
  # Endogenous variables
  vars_endogen <- unique(unlist(lapply(data, function(x) {dimnames(x[["endogen"]])[[2]]})))
  if (!is.null(variables)) {
    if (length(which(variables %in% vars_endogen)) == 0) {
      stop("Non of the variables specified in 'variables' is contained in the data.")
    }
    variables <- variables[which(variables %in% vars_endogen)]
  }
  
  data <- lapply(data, .diff_func, variables, multi)
  
  return(data)
}

.diff_func <- function(x, variables, multi){
  
  result <- x
  x <- result[["endogen"]]
  tsp_all <- stats::tsp(x)
  tsp_all[1] <- tsp_all[1] + 1 / tsp_all[3]
  if (is.null(multi)) {
    multi <- 1
  }
  if (is.null(variables)){
    for (i in dimnames(x)[[2]]){
      x[, i] <- c(NA, diff(x[, i])) * multi
    }
  } else {
    for (i in variables){
      if (is.element(i, dimnames(x)[[2]])){
        x[, i] <- c(NA, diff(x[, i])) * multi
      }
    } 
  }
  x <- x[-1, ]
  x <- stats::ts(x, start = tsp_all[1], frequency = tsp_all[3])
  
  result[["endogen"]] <- x
  return(result)
}