#' Plotting Forecasts of BVAR Models
#' 
#' A plot function for objects of class "bvarprd".
#' 
#' @param x an object of class "bvarprd", usually, a result of a call to \code{\link{predict.bgvar}}.
#' @param variable a character vector containing a country name as its first and a variable name
#' as its second element. If \code{NULL} (default) all series are plotted.
#' @param ... further graphical parameters.
#' 
#' @export
plot.bgvarprd <- function(x, variable = NULL, ...) {
  
  tt <- nrow(x[["y"]])
  if (is.null(variable)) {
    pos <- 1:ncol(x[["y"]])
  } else {
    pos <- which(x[["index"]][, "country"] == variable[1] & x$index[, "variable"] == variable[2])
  }
  for (i in pos) {
    n_ahead <- nrow(x[["fcst"]][[i]])
    temp <- matrix(NA, tt + n_ahead, 4)
    temp[1:tt, 1] <- x[["y"]][, i]
    temp[tt, 2:4] <- x[["y"]][tt, i]
    temp[(tt + 1):(tt + n_ahead), 2:4] <- x[["fcst"]][[i]]
    tsp_temp <- stats::tsp(x[["fcst"]][[i]])
    temp <- stats::ts(temp, end = tsp_temp[2], frequency = tsp_temp[3])
    rm(tsp_temp)
    title_temp <- paste(x[["index"]][i, "country"], ": ", x[["index"]][i, "variable"], sep = "")
    stats::plot.ts(temp, plot.type = "single", lty = c(1, 2, 1, 2), main = title_temp, ylab = "")
  }
  
}