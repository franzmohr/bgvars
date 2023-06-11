#' Posterior simulation
#' 
#' Estimates the country models of a Bayesian GVAR model.
#' 
#' @param object a list of data and model specifications for country-specific
#' models, which should be passed on to function \code{FUN}. Usually, the
#' output of a call to \code{\link{create_models}} in combination with \code{\link{add_priors}}.
#' @param FUN the function to be applied to each country model in argument \code{object}.
#' If \code{NULL} (default), the internal functions \code{\link{bvecxpost}} is used.
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' @param ctry character vector specifying for which countries posterior distributions
#' should be drawn. If \code{NULL} (default), draws are generated for all country models.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class \code{bgvecest}, which contains a list of data,
#' model specifications, priors, coefficient draws and information
#' criteria for each estimated country model.
#' 
#' @export
draw_posterior.gvecsubmodels <- function(object, ..., FUN = NULL, mc.cores = NULL, ctry = NULL){
  
  # If 'ctry' is specified, reduce list to relevant elements
  if (!is.null(ctry)) {
    pos <- which(names(object) %in% ctry)
    temp <- list()
    for (i in 1:length(pos)) {
      temp[[i]] <- object[[pos[i]]]
      names(temp)[i] <- names(object)[pos[i]]
    }
    rm(object)
    object <- temp
    rm(temp)
    names_obj <- names(object)
  }
  
  names_obj <- names(object)
  names_temp <- names_obj
  if (length(unique(names_temp)) != length(names_temp)) {
    for (i in unique(names_temp)) {
      pos_temp <- which(names_temp == i)
      temp <- names_temp[pos_temp]
      id <- paste0("0000", 1:length(temp))
      id <- substring(id, nchar(id) - 3, nchar(id))
      temp <- paste0(id, "-", temp)
      names_temp[pos_temp] <- temp
    }
  }
  names(object) <- names_temp
  
  # Print estimation information
  cat("Estimating submodels...\n")
  if (is.null(mc.cores)) {
    object <- lapply(object, .posterior_gvecsubmodels, use = FUN)
  } else {
    object <- parallel::mclapply(object, .posterior_gvecsubmodels, use = FUN,
                                 mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  
  names(object) <- names_obj
  
  if (is.null(ctry)) {
    class(object) <- append("bgvecest", class(object))
  }
  
  return(object)
}

# Helper function to implement try() functionality
.posterior_gvecsubmodels <- function(object, use) {
  
  if (is.null(use)) {
    object <- try(bvecxpost(object))
  } else {
    # Apply own function
    object <- try(use(object))
  }
  
  # Produce something if estimation fails
  if (inherits(object, "try-error")) {
    object <- c(object, list(error = TRUE))
  }
  
  return(object)
}