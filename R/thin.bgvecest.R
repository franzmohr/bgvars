#' Thinning Posterior Draws
#' 
#' Thins the MCMC posterior draws in an object of class \code{"bgvecest"}.
#' 
#' @param x an object of class \code{"bgvecest"}.
#' @param thin an integer specifying the thinning interval between successive values of posterior draws.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class \code{"bgvecest"}.
#' 
#' @export
thin.bgvecest <- function(x, thin = 10, ...) {
  
  vars <- c("a0", "a0_lambda", "a0_sigma",
            "alpha", "alpha_lambda", "alpha_sigma",
            "beta_domestic", "beta_domestic_lambda", "beta_domestic_sigma",
            "beta_foreign", "beta_foreign_lambda", "beta_foreign_sigma",
            "beta_global", "beta_global_lambda", "beta_global_sigma",
            "beta_deterministic", "beta_deterministic_lambda", "beta_deterministic_sigma",
            "gamma_domestic", "gamma_domestic_lambda", "gamma_domestic_sigma",
            "gamma_foreign", "gamma_foreign_lambda", "gamma_foreign_sigma",
            "gamma_global", "gamma_global_lambda", "gamma_global_sigma",
            "gamma_deterministic", "gamma_deterministic_lambda", "gamma_deterministic_sigma",
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