#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class 'submodelestlist'.
#' 
#' @param x an object of class 'submodelestlist'.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'submodelestlist'.
#' 
#' @export
thin.submodelestlist <- function(x, thin = 10, ...) {
  
  for (i in 1:length(x)) {
    
    if (!is.null(x[[i]][["error"]])) {
      if (x[[i]][["error"]]) {
        next
      }
    }
    
    x[[i]] <- thin(x[[i]], thin = thin, ...)
  }
  
  return(x)
}
