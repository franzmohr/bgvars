#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bgvarest"}.
#' 
#' @param x an object of class \code{"bgvarest"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class \code{"bgvarest"}.
#' 
#' @export
thin.bgvarest <- function(x, thin = 10, ...) {
  
  vars <- c("a0", "a0_lambda", "a0_sigma",
            "domestic", "domestic_lambda", "domestic_sigma",
            "foreign", "foreign_lambda", "foreign_sigma",
            "global", "global_lambda", "global_sigma",
            "deterministic", "deterministic_lambda", "deterministic_sigma",
            "sigma", "sigma_lambda")
  
  draws <- NA
  for (j in 1:length(x)) {
    
    if (!is.null(x[[j]][["error"]])) {
      if (x[[j]][["error"]]) {
        next
      }
    }
    
    for (i in vars) {
      if (is.na(draws)) {
        if (!is.null(x[[j]][["posteriors"]][[i]])) {
          if (is.list(x[[j]][["posteriors"]][[i]])) {
            draws <- nrow(x[[j]][["posteriors"]][[i]][[1]])
          } else {
            draws <- nrow(x[[j]][["posteriors"]][[i]]) 
          }
        }   
      }
    } 
  }
  
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  
  for (j in 1:length(x)) {
    
    if (!is.null(x[[j]][["error"]])) {
      if (x[[j]][["error"]]) {
        next
      }
    }
    
    for (i in vars) {
      if (!is.null(x[[j]][["posteriors"]][[i]])) {
        if (is.list(x[[j]][["posteriors"]][[i]])) {
          for (k in 1:length(x[[j]][["posteriors"]][[i]])) {
            x[[j]][["posteriors"]][[i]][[k]] <- coda::mcmc(as.matrix(x[[j]][["posteriors"]][[i]][[k]][pos_thin,]), start = start, end = end, thin = thin) 
          }
         } else {
          x[[j]][["posteriors"]][[i]] <- coda::mcmc(as.matrix(x[[j]][["posteriors"]][[i]][pos_thin,]), start = start, end = end, thin = thin)  
        }
      }
    } 
  }
  
  return(x)
}