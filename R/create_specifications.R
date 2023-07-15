#' Create Country Model Specifications
#' 
#' Produces a list of model specifications for each country in a GVAR model.
#' 
#' @param country_data a named list of time-series objects of country-specific data.
#' @param global_data a named time-series object of global data.
#' @param countries a character vector of country names, which should enter the GVAR model.
#' The names must correspond to the names of the elements in \code{country_data}. If
#' \code{NULL} (default), all countries are included.
#' @param domestic a named list of specificatios for domestic variables. See also 'Details'.
#' @param foreign a named list of specificatios for foreign variables. See also 'Details'.
#' @param global a named list of specificatios for global variables. See also 'Details'.
#' @param deterministic a named list of specifications for deterministic terms. See alos 'Details'.
#' @param type a character specifying if the country models are estimated as vector
#' autoregressive "VAR" (default) or vector error correction "VEC" models.
#' @param r an integer or vector for the cointegration rank of VEC models. See also 'Details'.
#' @param structural a character vector of country names, for which structural models
#' should be estimated. If \code{NULL} (default), no country-model will be structural.
#' @param tvp logical indicating whether the coefficients of the model are time varying (default is \code{FALSE}).
#' @param sv logical indicating whether the error variances are time varying and should be estimated using a stochastic volatility algorithm (default is \code{FALSE}).
#' @param iterations an integer of MCMC draws excluding burn-in draws (defaults
#' to 50000).
#' @param burnin an integer of MCMC draws used to initialize the sampler
#' (defaults to 5000). These draws do not enter the computation of posterior
#' moments, forecasts etc.
#' 
#' @details
#' 
#' For arguments \code{domestic}, \code{foreign} or \code{global} a list should be provided, which
#' contains a vector of the names of used variables in element \code{"variables"}, i.e., \code{variables = c("var1", "var2", ...)},
#' where \code{"var1"} is the name of a column in object \code{"country_data"} or \code{"global_data"}.
#' Similarly, the lag order of each variable group can be changed by adding an element \code{"lags"}
#' with an integer or integer vector, i.e. \code{lags = 1}. If a vector is provided,
#' \code{\link{create_models}} will produce a distinct model for all
#' possible specifications.
#' 
#' If a list is provided as argument \code{deterministic}, it can contain the following elements:
#' \describe{
#'  \item{\code{const}}{either logical for VAR models or a character specifying whether a constant term
#'  enters the error correction term (\code{"restricted"}) or the non-cointegration term (\code{"unrestricted"})
#'  of a VEC model. If `NULL` (default) no constant term will be added.}
#'  \item{\code{trend}}{either logical for VAR models or a character specifying whether a trend term
#'  enters the error correction term (\code{"restricted"}) or the non-cointegration term (\code{"unrestricted"})
#'  of a VEC model. If `NULL` (default) no trend term will be added.}
#'  \item{\code{seasonal}}{either logical for VAR models or a character specifying whether seasonal dummies
#'  enter the error correction term (\code{"restricted"}) or the non-cointegreation term (\code{"unrestricted"}).
#'  If `NULL` (default) no seasonal terms will be added. The amount of dummy variables depends
#'  on the frequency of the time series data in each list element.}
#' }
#' 
#' The argument \code{rank} is only used when \code{type = "VEC"}. Otherwise, the rank is set to zero.
#' If a vector is provided, \code{\link{create_models}} will produce a distinct model for all possible specifications.
#' 
#' @return A list of model specifications for each country.
#' 
#' @examples
#' # Load data
#' data("gvar2019")
#' 
#' # Create regions
#' temp <- create_regions(country_data = gvar2019$country_data,
#'              weight_data = gvar2019$weight_data,
#'              region_weights = gvar2019$region_weights,
#'              regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")),
#'              period = 3)
#' 
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' global_data = gvar2019$global_data
#' 
#' # Difference series to make them stationary
#' country_data <- diff_variables(country_data, variables = c("y", "Dp", "r"), multi = 100)
#' global_data <- diff_variables(global_data, multi = 100)
#' 
#' # Create time varying weights
#' weight_data <- create_weights(weight_data, period = 3, country_data = country_data)
#' 
#' # Generate specifications
#' model_specs <- create_specifications(
#'                  country_data = country_data,
#'                  global_data = global_data,
#'                  countries = c("US", "JP", "CA", "NO", "GB", "EA"), 
#'                  domestic = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  foreign = list(variables = c("y", "Dp", "r"), lags = 1),
#'                  global = list(variables = c("poil"), lags = 1),
#'                  deterministic = list(const = TRUE, trend = FALSE, seasonal = FALSE),
#'                  iterations = 10,
#'                  burnin = 10)
#' # Note that the number of iterations and burnin draws should be much higher!
#' 
#' @export
create_specifications <- function(country_data, global_data = NULL,
                                  countries = NULL,
                                  domestic,
                                  foreign,
                                  global = NULL,
                                  deterministic = NULL,
                                  type = "VAR",
                                  r = NULL,
                                  structural = NULL,
                                  tvp = FALSE,
                                  sv = FALSE,
                                  iterations = 50000,
                                  burnin = 5000){
  
  if (all(!c("VAR", "VEC") %in% type)) {"Argument 'type' must be either 'VAR' or 'VEC'."}
  
  # Basic data checks ----
  # Check if country data is a list
  if (!"list" %in% class(country_data)) {
    stop("Country data must be of class 'list'.")
  }
  # Check if data in country_data are of class ts
  if(any(!unlist(lapply(country_data, function(x) {"ts" %in% class(x)})))) {
    stop("Elements in 'country_data' must be of class 'ts'.")
  }
  # Check if global data is of class ts
  if (!is.null(global_data)) {
    if(!"ts" %in% class(global_data)) {
      stop("'global_data' must be of class 'ts'.")
    }
  }
  
  # Check if data are named
  if(is.null(names(country_data))){
    stop("'country_data' must be a named list. Please provide a name for each country.")
  }
  if(!is.null(global_data) & is.null(dimnames(global_data)[[2]])){
    stop("'global_data' must be named. Please, provide a name for each global variable.")
  }
  
  # Countries
  if (!is.null(countries)) {
    if (length(countries) < 2) {
      stop("Please, specifiy at least two countries.")
    }
  }
  
  use_structural <- FALSE
  if (!is.null(structural)) {
    
    if (!"character" %in% class(structural)) {
      stop("Argument 'structural' must be of class character.")
    }
    
    if (!is.null(countries)) {
      if (!all(structural %in% countries)) {
        stop("Not all countries in argument 'structural' are available.")
      }
    }
    use_structural <- TRUE
  }
  
  if (iterations <= 0) {
    stop("Argument 'iterations' must be a positive integer.")
  }
  if (burnin <= 0) {
    stop("Argument 'burnin' must be a positive integer.")
  }
  
  # Check deterministics ----
  if (!is.null(deterministic)) {
    
    if (is.null(names(deterministic))) {
      stop("Argument 'deterministic' must be a named list.")
    }
    
    if (!all(names(deterministic) %in% c("const", "trend", "seasonal"))) {
      stop("Invalid specification of argument 'deterministic'. Please check the names of the list entries.")
    }
    
    if (type == "VAR") {
      
      if ("const" %in% names(deterministic)) {
        if (!"logical" %in% class(deterministic[["const"]])) {
          stop("List element 'const' must be of class 'logical' for VAR models.")
        }
      }
      
      if ("trend" %in% names(deterministic)) {
        if (!"logical" %in% class(deterministic[["trend"]])) {
          stop("List element 'trend' must be of class 'logical' for VAR models.")
        }
      }
      
      if ("seasonal" %in% names(deterministic)) {
        if (!"logical" %in% class(deterministic[["seasonal"]])) {
          stop("List element 'seasonal' must be of class 'logical' for VAR models.")
        }
      }
    }
    
    if (type == "VEC") {
      if (!is.null(deterministic[["const"]]) & "logical" %in% class(deterministic[["const"]])) {
        stop("List element 'const' must not be logical for VEC models.")
      }
      if (!is.null(deterministic[["trend"]]) & "logical" %in% class(deterministic[["trend"]])) {
        stop("List element 'trend' must not be logical for VEC models.")
      }
      if (!is.null(deterministic[["seasonal"]]) & "logical" %in% class(deterministic[["seasonal"]])) {
        stop("List element 'seasonal' must not be logical for VEC models.")
      }
    }
  }
  
  variables <- NULL
  p_domestic <- 1
  if (!is.null(domestic)) {
    if (!is.null(domestic$variables)) {
      variables <- c(variables, domestic$variables)
    }
    if (!is.null(domestic$lags)) {
      if (type == "VEC" & 0 %in% domestic$lags) {
        stop("The lag order of domestc variables must be at least 1 for VAR models.")
      }
      p_domestic <- domestic$lags
    }
  }
  
  p_foreign <- 1
  if (!is.null(foreign)) {
    if (!is.null(foreign$variables)) {
      variables <- c(variables, foreign$variables)
    }
    if (!is.null(foreign$lags)) {
      if (type == "VEC" & 0 %in% foreign$lags) {
        stop("The lag order of foreign variables must be at least 1 for VEC models.")
      }
      p_foreign <- foreign$lags
    }
  } else {
    variables <- unique(c(variables, unlist(lapply(country_data, function(x) {dimnames(x)[[2]]}))))
  }
  
  s <- 1
  if (!is.null(global)) {
    if (!is.null(global$variables)) {
      variables <- c(variables, global$variables)
    }
    if (!is.null(global$lags)) {
      if (type == "VEC" & 0 %in% global$lags) {
        stop("The lag order of global variables must be at least 1 for VEC models.")
      }
      s <- global$lags
    }
  }
  
  # Prepare specifications ----
  if (any(p_foreign < 1)) {
    stop("Lag of foreign star variables in VARX model must be at least 1.")
  }

  s.vars <- unique(unlist(lapply(country_data, function(x) {return(dimnames(x)[[2]])})))
  if (!is.null(variables)){
    variables <- unique(variables)
    s.vars <- s.vars[which(is.element(s.vars, variables))]
  }
  
  use_global <- FALSE
  if (is.null(global)){
    g.vars <- NA
  } else {
    if (is.null(global_data)) {
      stop("Argument 'global' was specified, but no data was provided in 'global_data'.")
    }
    g.vars <- dimnames(global_data)[[2]]
    if (!is.null(global$variables)) {
      g.used <- which(is.element(g.vars, global$variables)) 
    } else {
      g.used <- 1:length(g.vars)
    }
    if (length(g.used) > 0){
      g.vars <- g.vars[g.used]
      use_global <- TRUE
    }
  }
  global <- FALSE
  if (use_global) {
    global <- TRUE
  }
  
  # Deterministic terms ----
  use_det <- FALSE
  if (type == "VAR" & !is.null(deterministic)) {
    use_det <- TRUE
    det_var <- NULL
    if ("const" %in% names(deterministic)) { # Check if specification was provided
      if (deterministic$const) { # Check if TRUE
        det_var <- c(det_var, "const") 
      }
    }
    if ("trend" %in% names(deterministic)) { # Check if specification was provided
      if (deterministic$trend) {# Check if TRUE
        det_var <- c(det_var, "trend") 
      }
    }
    if ("seasonal" %in% names(deterministic)) { # Check if specification was provided
      if (deterministic$seasonal) {# Check if TRUE
        det_var <- c(det_var, "seasonal")
      }
    }
  }
  if (type == "VEC" & !is.null(deterministic)) {
    use_det <- TRUE
    det_res <- NULL
    det_unres <- NULL
    if ("const" %in% names(deterministic)) { # Check if specification was provided
      if (deterministic$const == "restricted") { # Add restricted dets
        det_res <- c(det_res, "const") 
      }
      if (deterministic$const == "unrestricted") { # Add unrestricted dets
        det_unres <- c(det_unres, "const") 
      }
    }
    if ("trend" %in% names(deterministic)) { # Check if specification was provided
      if (deterministic$trend == "restricted") { # Add restricted dets
        det_res <- c(det_res, "trend") 
      }
      if (deterministic$trend == "unrestricted") { # Add unrestricted dets
        det_unres <- c(det_unres, "trend") 
      }
    }
    if ("seasonal" %in% names(deterministic)) { # Check if specification was provided
      if (deterministic$seasonal == "restricted") { # Add restricted dets
        det_res <- c(det_res, "seasonal") 
      }
      if (deterministic$seasonal == "unrestricted") { # Add unrestricted dets
        det_unres <- c(det_unres, "seasonal") 
      }
    }
  }
  
  # Create list of specifications ----
  specs <- list()
  specs.names <- NULL
  for (i in 1:length(country_data)){
    d.vars <- dimnames(country_data[[i]])[[2]]
    if (!is.null(domestic$variables)){
      d.vars <- d.vars[which(is.element(d.vars, domestic$variables))]
    }
    s.vars_i <- s.vars
    if (!is.null(foreign$variables)) {
      s.vars_i <- foreign$variables
    }
    specs <- c(specs, list(list(type = type,
                                domestic = list("variables" = d.vars, "lags" = p_domestic),
                                foreign = list("variables" = s.vars_i, "lags" = p_foreign))))
    if (global) {
      specs[[i]]$global <- list("variables" = g.vars, "lags" = s)
    }
    
    if (use_det) {
      if (type == "VAR") {
        specs[[i]]$deterministic <- det_var
      }
      if (type == "VEC") {
        if (!is.null(det_res)) {
          specs[[i]]$deterministic$restricted <- det_res
        }
        if (!is.null(det_unres)) {
          specs[[i]]$deterministic$unrestricted <- det_unres
        }
      }
    }
    
    if (type == "VEC") {
      if (is.null(r)) {
        specs[[i]][["rank"]] = 0:length(d.vars)
      } else {
        specs[[i]][["rank"]] = r
      }
    }
    
    specs[[i]][["structural"]] <- FALSE  
    if (use_structural) {
      if (names(country_data)[i] %in% structural) {
        specs[[i]][["structural"]] <- TRUE
      }
    }
    
    specs[[i]][["tvp"]] <- tvp
    specs[[i]][["sv"]] <- sv
    specs[[i]][["iterations"]] <- iterations
    specs[[i]][["burnin"]] <- burnin
    specs.names <- c(specs.names, names(country_data)[i])
  }
  names(specs) <- specs.names
  for (i in 1:length(specs)){
    specs[[i]]$domestic$variables <- specs[[i]]$domestic$variables[is.element(specs[[i]]$domestic$variables, dimnames(country_data[[i]])[[2]])]
  }
  if (!is.null(countries)){
    if (any(!is.element(countries, names(country_data)))) {
      stop(paste("There are no country data available for ", paste(countries[which(!is.element(countries, names(country_data)))], collapse = ", "), ".", sep = ""))
    }
    result <- NULL
    names.result <- NULL
    for (i in countries){
      result <- c(result, list(specs[[i]]))
      names.result <- c(names.result, i)
    }
    names(result) <- names.result
  } else {
    result <- specs
  }
  return(result)
}