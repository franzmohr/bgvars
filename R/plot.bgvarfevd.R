#' Plotting Generalised Forecast Error Variance Decompositions of Bayesian Global Vector Autoregressions
#' 
#' A plot function for objects of class \code{"bgvarfevd"}.
#' 
#' @param x an object of class \code{"bgvarfevd"}, usually, a result of a call to \code{\link{gfevd}}.
#' @param top_n integer specifying the amount of the most explaining variables that should be displayed
#' @param group either \code{"variable"} (default) for variable-wise or \code{"country"} for
#' country-wise decompositions.
#' @param ... further graphical parameters.
#' 
#' @export
#' @rdname gfevd
plot.bgvarfevd <- function(x, top_n = 5L, group = "variable", ...) {
  x_temp <- t(apply(x, 1, function(x) {x / sum(x)}))
  if (group == "country") {
    c_names <- strsplit(dimnames(x_temp)[[2]], "_")
    c_names <- unlist(lapply(c_names, function(x) {x[1]}))
    countries <- unique(c_names)
    temp <- NULL
    for (i in countries) {
      temp <- cbind(temp, rowSums(x[, which(i == c_names)]))
    }
    dimnames(temp) <- list(NULL, countries)
    x <- stats::ts(temp, start = 0, frequency = 1)
  }
  n_ctr <- ncol(x)
  x_temp <- colSums(x)
  top <- order(-x_temp)[1:top_n] 

  row <- x[, -top]
  if (n_ctr - length(top) > 1) {
    row <- rowSums(row)
  } else {
    row <- matrix(row)
  }
  result <- cbind(x[, top], row)
  
  if (group == "variable") {
    legend_names <- strsplit(dimnames(x)[[2]][top], "_")
    for (i in 1:length(legend_names)) {
      legend_names[[i]] <- paste(legend_names[[i]][1], ": ", legend_names[[i]][2], sep = "")
    } 
    legend_names <- unlist(legend_names)
  }
  if (group == "country") {
    legend_names <- dimnames(x)[[2]][top]
  }
  legend_names <- c(legend_names[1:top_n], "RoW")
  
  par_orig <- graphics::par("mar")
  graphics::par(mar = c(5.1, 4.1, 4.1, 8.1))
  graphics::barplot(t(result), ylab = "Percentage", xlab = "Period", names.arg = stats::time(result))
  graphics::par(mar = par_orig)
  graphics::legend("left", legend = legend_names, xpd = FALSE, fill = grDevices::gray.colors(top_n + 1), inset = 1)
}
