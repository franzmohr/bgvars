#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bgvar"}.
#' 
#' @param x an object of class \code{"bgvar"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class \code{"bgvar"}.
#' 
#' @export
thin.bgvar <- function(x, thin = 10, ...) {
  
  draws <- NA
  vars <- c("a0", "a", "b", "c", "sigma")
  for (i in vars) {
    if (is.na(draws)) {
      if (!is.null(x[[i]])) {
        if (is.list(x[[i]])) {
          draws <- nrow(x[[i]][[1]])
        } else {
          draws <- nrow(x[[i]]) 
        }
      }   
    }
  }
  
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  for (i in vars) {
    if (!is.null(x[[i]])) {
      x[[i]] <- coda::mcmc(as.matrix(x[[i]][pos_thin,]), start = start, end = end, thin = thin)  
    }
  }
  
  return(x)
}