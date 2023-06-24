#' Plotting Impulse Responses of Bayesian GVAR Submodels
#' 
#' A plot function for objects of class \code{"bgvarestirf"}.
#' 
#' @param x an object of class \code{"bgvarestirf"}, usually, a result of a call to \code{\link{irf.bgvarest}}.
#' @param ... further graphical parameters.
#' 
#' @export 
#' @rdname irf.bgvarest
plot.bgvarestirf <- function(x, ...) {
  for (i in 1:length(x)) {
    
    # Prepare title and axis text
    main <- attr(x[[i]], "ctry")
    if (!is.null(attr(x[[i]], "rank"))) {
      main <- paste0(main, " (r = ", attr(x[[i]], "rank"),")")
    }
    ylab <- paste0(attr(x[[i]], "impulse"), " -> ", attr(x[[i]], "response"))
    
    # Plot
    plot(x[[i]], main = main, ylab = ylab, xlab = NULL, ...)
    
    rm(list = c("main", "ylab")) # Housekeeping
    
  }
}


