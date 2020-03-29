#' Fitter Function for GVAR Models
#' 
#' Estimates the country models of a Bayesian GVAR model.
#' 
#' @param object a list of data and model specifications for country-specific
#' models and information on the global VAR model as produced by the
#' \code{country_models} function.
#' @param iterations an integer of MCMC draws including burn-in (defaults
#' to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 5000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' @param thin an integer specifying the thinning factor for the MCMC output.
#' Defaults to 10, which means that the forecast sequences contain only every
#' tenth draw of the original sequence. Set \code{thin = 1} to obtain the
#' full MCMC sequence.
#' @param mc.cores the number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. The option is initialized from
#' environment variable MC_CORES if set. Must be at least one, and
#' parallelization requires at least two cores.
#' 
#' @return An object of class \code{bgvarest} , which contains a list of data,
#' model specifications, priors, (thinned) coefficient draws and information
#' criteria for each estimated country model.
#' 
#' @examples 
#' data("gvar2016")
#' 
#' country_data <- gvar2016$country_data
#' global_data <- gvar2016$global_data
#' region_weights <- gvar2016$region_weights
#' weight_data <- gvar2016$weight_data
#' 
#' # Take first difference of all variables y and Dp
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"))
#' 
#' # Generate EA area region with 2 year, rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 2,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 2 year rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "r"),
#'                                      countries = c("EA", "US", "JP", "CA", "GB"),
#'                                      p_domestic = 1,
#'                                      p_foreign = 1,
#'                                      type = "VAR")
#' 
#' # Create estimable objects
#' object <- create_models(country_data = country_data,
#'                         gvar_weights = gvar_weights,
#'                         model_specs = model_specs)
#' 
#' # Add priors
#' object <- add_priors(object)
#' 
#' # Estimate GVAR model
#' gvar_est <- estimate_gvar(object, iterations = 50, burnin = 10, thin = 4)
#' # Note that the number of iterations and burnin should be much higher.
#' 
#' @export
estimate_gvar <- function(object, iterations = 50000, burnin = 5000, thin = 10, mc.cores = NULL){
  if (iterations <= 0) {
    stop("Argument 'iterations' must be a positive integer.")
  }
  if (burnin <= 0) {
    stop("Argument 'burnin' must be a positive integer.")
  }
  if (thin < 1) {
    stop("Argument 'thin' must be at least 1.")
  }
  
  post_burnin <- iterations - burnin
  
  if (post_burnin <= 0) {
    stop("Argument 'iterations' must be larger than argument 'burnin'.")
  }
  
  if (post_burnin < thin) {
    stop("Not enough post-burn-in draws for the specified thinning factor. Consider increasing the value of argument 'iterations'.")
  }
  
  # Check if specified number of draws for variable selection and argument thin are valid
  vs <- unlist(lapply(object, function(x){!is.null(x$variable_selection)}))
  if (any(vs)) {
    threshold <- unlist(lapply(object, function(x){!is.null(x$variable_selection$threshold)}))
    if (any(threshold)) {
      for (i in which(threshold)) {
        if (!is.null(object[[i]]$variable_selection$draws)) {
          if (object[[i]]$variable_selection$draws > post_burnin) {
            stop("Number of variable selection draws may not be larger than the number of post-burn-in draws.")
          }
        }
      }
    }
  }
  
  # Print estimation information
  cat("Estimating country models.\n")
  if (is.null(mc.cores)) {
    object <- lapply(object, .estimate_country_model, iterations = iterations, burnin = burnin, thin = thin)
  } else {
    object <- parallel::mclapply(object, .estimate_country_model,
                                 iterations = iterations, burnin = burnin, thin = thin,
                                 mc.cores = mc.cores, mc.preschedule = FALSE)
  }
  class(object) <- append("bgvarest", class(object))
  return(object)
}

.estimate_country_model <- function(object, iterations, burnin, thin) {
  type <- object$model$type
  if (type == "VAR"){
    object <- try(.bvarx(object, iterations = iterations, burnin = burnin, thin = thin))
  }
  if (type == "VEC") {
    object <- try(.bvecx_rr(object, iterations = iterations, burnin = burnin, thin = thin))
  }
  if (inherits(object, "try-error")) {
    object <- c(object, list(coefficients = NULL))
  }
  return(object)
}