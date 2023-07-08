#' Plotting Submodel Test Statistics
#' 
#' A plot function for objects of class \code{"ctryvartest"} for visual inspection
#' of posterior draws.
#' 
#' @param x an object of class \code{"ctryvarest"}, usually, a result of a call
#' to \code{\link{draw_posterior}}.
#' @param ctry character vector (optional). Names of the countries, for which test
#' statistics should be visualized.
#' @param ... further graphical parameters.
#' 
#' @export
plot.teststats.bgvarest <- function(x, ctry = NULL, ...) {
  
  if (is.null(ctry)) {
    pos <- 1:length(x)
  } else {
    pos <- which(names(x) %in% ctry)
  }
  
  mgp_orig <- graphics::par("mgp")
  graphics::par(mgp = c(3, 2.5, 0))
  
  for (i in pos) {
    ctry_name <- names(x)[i]
    ctry_temp <- c()
    x_names <- c()
    for (j in 1:nrow(x[[i]][["teststats"]])) {
      
      if (is.na(x[[i]][["teststats"]][j, "LL"])) {
        next
      }
      
      ctry_temp[[j]] <- colSums(x[[i]][[2]][[j]])
      
      temp <- c()
      if ("p_domestic" %in% names(x[[i]][["teststats"]])) {
        temp <- append(temp, paste0("p = ", x[[i]][["teststats"]][j, "p_domestic"]))
      }
      if ("p_foreign" %in% names(x[[i]][["teststats"]])) {
        temp <- append(temp, paste0("p* = ", x[[i]][["teststats"]][j, "p_foreign"]))
      }
      if ("s" %in% names(x[[i]][["teststats"]])) {
        temp <- append(temp, paste0("s = ", x[[i]][["teststats"]][j, "s"]))
      }
      if ("r" %in% names(x[[i]][["teststats"]])) {
        temp <- append(temp, paste0("r = ", x[[i]][["teststats"]][j, "r"]))
      }
      temp <- paste0(temp, collapse = "\n")
      x_names <- append(x_names, temp)
    }
    
    names(ctry_temp) <- x_names
    graphics::boxplot(ctry_temp, main = ctry_name, ...)
  }
  
  graphics::par(mgp = mgp_orig)
  
}


