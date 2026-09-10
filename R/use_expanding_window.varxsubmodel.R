#' Expanding Window Estimation
#' 
#' Creates objects for expanding window posterior simulation.
#' 
#' @param object a list of model specifications. Usually, the output of a call to 
#' a model creation function in combination with \code{\link[bvartools]{add_priors}} and
#' \code{\link[bvartools]{add_initial_values}}.
#' @param start the start period of the prediction of the first iteration of the
#' expanding window approach.
#' @param ... arguments passed forward to method.
#' 
#' @export
use_expanding_window.varxsubmodel <- function(object, start, ...) {

  k <- object[["model"]][["k_domestic"]]
  y <- object[["data"]][["y"]]
  tsp_y <- stats::tsp(y)
  time_y <- stats::time(y)
  test <- stats::window(y, start = start)
  time_train <- time_y[time_y < min(stats::time(test))]
  nobs_train_min <- length(time_train)
  nobs_train_max <- length(time_y)
  pos_end <- nobs_train_min:nobs_train_max
  
  result <- list()
  for (i in 1:length(pos_end)) {
    
    temp <- object
    
    # Trim data
    dims_y <- dimnames(temp[["data"]][["y"]])
    temp[["data"]][["y"]] <- stats::ts(as.matrix(temp[["data"]][["y"]][1:pos_end[i], ]),
                                       start = tsp_y[1], frequency = tsp_y[3], class = c("mts", "ts", "matrix"))
    dimnames(temp[["data"]][["y"]]) <- dims_y
    if (!is.null(temp[["data"]][["x"]])) {
      dims_x <- dimnames(temp[["data"]][["x"]])
      temp[["data"]][["x"]] <- stats::ts(as.matrix(temp[["data"]][["x"]][1:pos_end[i], ]),
                                         start = tsp_y[1], frequency = tsp_y[3], class = c("mts", "ts", "matrix"))
      dimnames(temp[["data"]][["x"]]) <- dims_x
    }
    if (!is.null(temp[["data"]][["z"]])) {
      temp[["data"]][["z"]] <- temp[["data"]][["z"]][1:(k * pos_end[i]), ]
    }
    
    result[[i]] <- temp
    rm(temp)
  }
  
  class(result) <- append("expwinmodellist", class(result))
  
  return(result)
}