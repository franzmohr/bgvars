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


#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class 'varxsubmodelest'.
#' 
#' @param x an object of class 'varxsubmodelest'.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'varxsubmodelest'.
#' 
#' @export
thin.varxsubmodelest <- function(x, thin = 10, ...) {
  
  vars <- c("sigma", "a")
  
  draws <- NA
  for (i in vars) {
    if (is.na(draws)) {
      if (!is.null(x[["posteriors"]][[i]])) {
        if (is.list(x[["posteriors"]][[i]][["coeffs"]])) {
          draws <- nrow(x[["posteriors"]][[i]][["coeffs"]][[1]])
        } else {
          draws <- nrow(x[["posteriors"]][[i]][["coeffs"]]) 
        }
      }   
    }
  } 
  
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  x[["model"]][["iterations"]] <- length(pos_thin)
  
  for (i in vars) {
    if (!is.null(x[["posteriors"]][[i]])) {
      if (is.list(x[["posteriors"]][[i]][["coeffs"]])) {
        for (k in 1:length(x[["posteriors"]][[i]][["coeffs"]])) {
          x[["posteriors"]][[i]][["coeffs"]][[k]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["coeffs"]][[k]][pos_thin,]), start = start, end = end, thin = thin) 
        }
      } else {
        x[["posteriors"]][[i]][["coeffs"]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["coeffs"]][pos_thin,]), start = start, end = end, thin = thin)  
      }
      if (!is.null(x[["posteriors"]][[i]][["lambda"]])) {
        x[["posteriors"]][[i]][["lambda"]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["lambda"]][pos_thin,]), start = start, end = end, thin = thin)  
      }
      if (!is.null(x[["posteriors"]][[i]][["sigma"]])) {
        x[["posteriors"]][[i]][["sigma"]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["sigma"]][pos_thin,]), start = start, end = end, thin = thin)  
      }
    }
  } 
  
  return(x)
}




#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class 'vecxsubmodelest'.
#' 
#' @param x an object of class 'vecxsubmodelest'.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'vecxsubmodelest'.
#' 
#' @export
thin.vecxsubmodelest <- function(x, thin = 10, ...) {
  
  vars <- c("sigma", "a", "beta")
  
  draws <- NA
  for (i in vars) {
    if (is.na(draws)) {
      if (!is.null(x[["posteriors"]][[i]])) {
        if (is.list(x[["posteriors"]][[i]][["coeffs"]])) {
          draws <- nrow(x[["posteriors"]][[i]][["coeffs"]][[1]])
        } else {
          draws <- nrow(x[["posteriors"]][[i]][["coeffs"]]) 
        }
      }   
    }
  } 
  
  pos_thin <- seq(from = thin, to = draws, by = thin)
  start <- pos_thin[1]
  end <- pos_thin[length(pos_thin)]
  x[["model"]][["iterations"]] <- length(pos_thin)
  
  for (i in vars) {
    if (!is.null(x[["posteriors"]][[i]])) {
      if (is.list(x[["posteriors"]][[i]][["coeffs"]])) {
        for (k in 1:length(x[["posteriors"]][[i]][["coeffs"]])) {
          x[["posteriors"]][[i]][["coeffs"]][[k]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["coeffs"]][[k]][pos_thin,]), start = start, end = end, thin = thin) 
        }
      } else {
        x[["posteriors"]][[i]][["coeffs"]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["coeffs"]][pos_thin,]), start = start, end = end, thin = thin)  
      }
      if (!is.null(x[["posteriors"]][[i]][["lambda"]])) {
        x[["posteriors"]][[i]][["lambda"]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["lambda"]][pos_thin,]), start = start, end = end, thin = thin)  
      }
      if (!is.null(x[["posteriors"]][[i]][["sigma"]])) {
        x[["posteriors"]][[i]][["sigma"]] <- coda::mcmc(as.matrix(x[["posteriors"]][[i]][["sigma"]][pos_thin,]), start = start, end = end, thin = thin)  
      }
    }
  }
  
  return(x)
}
