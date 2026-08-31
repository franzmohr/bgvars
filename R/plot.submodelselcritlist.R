#' Plotting Selection Criteria
#' 
#' A plot function for objects of class 'submodelselcritlist'.
#' 
#' @param x an object of class 'submodelselcritlist', usually, a result of a call
#' to \code{\link{selection_criteria}}.
#' @param criterion the selection criterion that should be plotted. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"} (default), \code{"HQ"}.
#' @param ... further graphical parameters passed on to \link[graphics]{boxplot}.
#' 
#' @export
plot.submodelselcritlist <- function(x, criterion = "BIC", ...) {
  
  submodels <- unique(names(x))
  
  for (i in submodels) {
    temp_list <- list()
    pos_model <- which(names(x) == i)
    for (j in 1:length(pos_model)) {
      temp_list[[j]] <- x[[pos_model[j]]]
    }
    class(temp_list) <- append("selcritlist", class(temp_list))
    plot(temp_list, criterion = criterion, main = i, ...)
  }
  
}
