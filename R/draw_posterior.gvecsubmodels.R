#' Posterior simulation
#' 
#' Estimates the country models of a Bayesian GVAR model.
#' 
#' @param object a list of data and model specifications for country-specific
#' models, which should be passed on to function \code{FUN}. Usually, the
#' output of a call to \code{\link{create_models}} in combination with \code{\link{add_priors}}.
#' @param FUN the function to be applied to each country model in argument \code{object}.
#' If \code{NULL} (default), the internal functions \code{\link{countryvarx}} is used for
#' VAR models and \code{\link{countryvecx}} for VEC models.
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
#' @examples 
#' 
#' # Load gvar2016 dataset
#' data("gvar2016")
#'
#' # Data objects
#' country_data <- gvar2016[["country_data"]]
#' global_data <- gvar2016[["global_data"]]
#' region_weights <- gvar2016[["region_weights"]]
#' weight_data <- gvar2016[["weight_data"]]
#' 
#' # Obtain weights
#' weight_data <- create_weights(weight_data = weight_data, period = 3,
#'                               country_data = country_data)
#' 
#' # Generate specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      countries = c("US", "JP", "CA", "NO", "GB"), 
#'                                      domestic = list(variables = c("y", "Dp", "r"), lags = 1),
#'                                      foreign = list(variables = c("y", "Dp", "r"), lags = 1),
#'                                      global = list(variables = c("poil"), lags = 1),
#'                                      deterministic = list(const = TRUE, trend = FALSE, seasonal = FALSE),
#'                                      iterations = 10,
#'                                      burnin = 10,
#'                                      thin = 1)
#'                                      
#' # Note that the number of iterations and burnin draws is usually much higher.
#' 
#' # Create list element for each model
#' country_models <- create_models(country_data = country_data,
#'                                 weight_data = weight_data,
#'                                 global_data = global_data,
#'                                 model_specs = model_specs)
#'                                 
#' # Add priors
#' country_models <- add_priors(country_models)
#' 
#' # Draw from posterior distribution
#' gvar_est <- draw_posterior(country_models, mc.cores = 1)
#' 
#' @export
draw_posterior.gvecsubmodels <- function(object, FUN = NULL, mc.cores = NULL, ctry = NULL){
  
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
    object <- c(object, list(coefficients = NULL))
  }
  
  return(object)
}