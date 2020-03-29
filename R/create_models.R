#' Create Country Models
#'
#' Produces a list of country models for each model specification that should be estimated.
#'
#' @param country_data a named list of time-series objects of country-specific data.
#' @param gvar_weights a matrix or an array of weights, usually, the output of a call to \code{\link{create_weights}}.
#' @param global_data a named time-series object of global data.
#' @param model_specs a list of model specifications as produced by \code{\link{create_specifications}}.
#' 
#' @return A list of country model input.
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
#' country_data <- diff_variables(country_data, variables = c("y", "Dp"))
#' 
#' # Generate weight matrices as 2 year, rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "poil"),
#'                                      p_domestic = 2,
#'                                      p_foreign = 2,
#'                                      type = "VAR")
#' 
#' # Create estimable objects
#' object <- create_models(country_data = country_data,
#'                         gvar_weights = gvar_weights,
#'                         global_data = global_data,
#'                         model_specs = model_specs)
#' 
#' @export
create_models <- function(country_data, gvar_weights, global_data = NULL,
                          model_specs = NULL){
  
  # Check if weights are named
  if (is.null(dimnames(gvar_weights)[[1]]) | is.null(dimnames(gvar_weights)[[2]])) {
    stop("Rows and columns of 'gvar_weights' must both be named.")
  }
  
  #### Prepare weight matrices ####
  countries <- names(model_specs)
  const_weight <- class(gvar_weights) == "matrix"
  if (const_weight) {
    temp_w <- matrix(0, length(countries), length(countries))
    dimnames(temp_w) <- list(countries, countries)
  } else {
    temp_w <- array(0, dim = c(length(countries), length(countries), dim(gvar_weights)[[3]]))
    dimnames(temp_w) <- list(countries, countries, dimnames(gvar_weights)[[3]])
  }
  
  # Recalculate weights
  temp_c <- NULL
  for (i in countries){
    temp_c <- c(temp_c, list(country_data[[i]]))
    if (const_weight){
      temp_w[i,] <- gvar_weights[i, countries] / (1 - sum(gvar_weights[i, !is.element(dimnames(gvar_weights)[[2]], countries)]))
    } else {
      for (j in 1:dim(gvar_weights)[3]){
        temp_w[i,,j] <- gvar_weights[i, countries, j] / (1 - sum(gvar_weights[i, !is.element(dimnames(gvar_weights)[[2]], countries), j])) 
      }
    }
  }
  country_data <- temp_c; rm(temp_c)
  names(country_data) <- countries
  gvar_weights <- temp_w; rm(temp_w)
  
  # Check if all country time series have the same length
  n_c <- length(countries)
  t_max <- c()
  for (i in 1:n_c){
    if (any(apply(country_data[[i]], 2, function(x){any(is.na(x))}))){
      vars <- dimnames(country_data[[i]])[[2]][which(apply(country_data[[i]], 2, function(x){any(is.na(x))}))]
      country <- names(country_data)[i]
      stop(paste("In country ", country, " the variable(s) ", paste(vars, sep = ", "), " contain NA values.", sep = ""))
      rm(vars)
    }
    t_max <- append(t_max, dim(stats::na.omit(country_data[[i]]))[1])
  }
  
  if (min(t_max) != max(t_max)) {
    stop("Numbers of observations differ across countries. For now, this package requires the same number for each country.")
  } else {
    t <- t_max[1]
    rm(t_max)
  }
  
  # Check if weight matrix has the same order of names in the rows and columns
  if (any(dimnames(gvar_weights)[[1]] != dimnames(gvar_weights)[[2]])){
    stop("The order of country names in the rows and columns of the weight matrix is not the same.")
  }
  
  # Check if weight matrix contains data on all required countries
  if (any(!is.element(names(country_data), dimnames(gvar_weights)[[1]]))){
    stop("The weight matrix does not contain data for at least one country in the country sample or is named differently.")
  }
  
  # If weight matrix is time varying, check if the number of periods is the same as in the country series
  if (!is.na(dim(gvar_weights)[3])){
    if (dim(gvar_weights)[3] > t) {
      warning("The weight array does not contain as much periods as the country sample. Trying to correct this.")
      dim_w <- as.numeric(dimnames(gvar_weights)[[3]])
      w_t <- stats::ts(dim_w, start = dim_w[1], frequency = stats::frequency(country_data[[1]]))
      gvar_weights <- gvar_weights[,,-which(!w_t %in% stats::time(country_data[[1]]))]
    }
  }
  
  #### Reduce object with global data to relevant series  ####
  global_vars <- any(unlist(lapply(model_specs, function(x) {return(!is.null(x$global))})))
  if (global_vars){
    if (is.null(global_data)) {
      stop("The model specifications contain global variables, but argument 'global_data' was not specified.")
    }
    global_vars_temp <- unique(unlist(lapply(model_specs, function(x) {return(x$global$variables)})))
    tsp_temp <- stats::tsp(global_data)
    global_data <- stats::as.ts(as.matrix(global_data[, global_vars_temp]))
    dimnames(global_data)[[2]] <- global_vars_temp
    stats::tsp(global_data) <- tsp_temp
    rm(tsp_temp)
  }
  
  if (!is.null(global_data) & !global_vars) {
    stop("Global data were provided in 'global_data', but they are not contained in the model specifications. Please adapt specifications, if you want to use global variables.")
  }
  rm(global_vars)
  
  #### Create global model data ####
  # Get list of all variables
  avail.dom.vars <- unique(unlist(lapply(country_data, function(x) {return(dimnames(x)[[2]])}))) # Domestic variables
  if (!is.null(global_data)){
    avail.global.vars <- dimnames(global_data)[[2]] # Global variables
  }
  
  # Generate final variable index
  index <- data.frame() # Index
  X <- c() # Final data set of domestic variables (incl. global endogenous)
  X_names <- c()
  for (i in countries){
    temp <- c()
    temp_names <- c()
    for (j in model_specs[[i]]$domestic$variables){
      if (is.element(j, dimnames(country_data[[i]])[[2]])){
        temp <- cbind(temp, country_data[[i]][, j])
      } else {
        temp <- cbind(temp, global_data[, j])
      }
      temp_names <- append(temp_names, j)
    }
    X <- cbind(X, temp)
    n_temp <- length(temp_names)
    X_names <- c(X_names, temp_names)
    index <- rbind(index, data.frame("Country" = rep(i, n_temp),
                                     "Variable" = temp_names))
  }
  X <- stats::na.omit(X)
  dimnames(X)[[2]] <- X_names
  index[, 1] <- as.character(index[, 1])
  index[, 2] <- as.character(index[, 2])
  
  #### Update modelled global variables ####
  variables <- unique(index[, 2]) # All endogenous variables
  endo.g.var <- NULL
  for (i in countries){
    # Check if there is an endogenous global variable in a country specification
    if (any(is.element(model_specs[[i]]$global$variables, variables))){
      foreign <- model_specs[[i]]$foreign$variables
      global <- model_specs[[i]]$global$variables
      # Which global variable is endogenous
      pos.endo.g <- which(is.element(global, variables))
      # Endo global vars
      endo.global <- global[pos.endo.g]
      # Add to vector of endo global vars
      endo.g.var <- c(endo.g.var, endo.global)
      # If a global variable is endogenous in another country model
      # add it to the foreign variable vector
      if (!any(endo.global %in% model_specs[[i]]$domestic$variables)) {
        foreign <- c(foreign, endo.global)
        foreign <- unique(foreign)
      }
      # Remove from global variables
      global <- global[-pos.endo.g]
      if (length(global) == 0){
        model_specs[[i]]$global <- NULL
      } else {
        model_specs[[i]]$global$variables <- global
      }
      model_specs[[i]]$foreign$variables <- foreign
    }
  }
  # Update global data
  endo.g.var <- unique(endo.g.var)
  if (length(endo.g.var) > 0){
    if (length(endo.g.var) == length(dimnames(global_data)[[2]])){
      global_data <- NULL
    } else {
      global_data <- global_data[, -which(is.element(endo.g.var, dimnames(global_data)[[2]]))]
    }
  }
  
  #### Trim data to same length ####
  global <- !is.null(global_data) # If a global variable is used
  if (global){
    k <- dim(X)[2]
    names_X <- dimnames(X)[[2]]
    names_g <- dimnames(global_data)[[2]]
    temp <- stats::na.omit(cbind(X, global_data))
    X <- temp[, 1:k]
    X_global <- stats::as.ts(as.matrix(temp[, -(1:k)]))
    stats::tsp(X_global) <- stats::tsp(X)
    dimnames(X)[[2]] <- names_X
    dimnames(X_global)[[2]] <- names_g
    rm(list = c("temp", "names_X", "names_g"))
  } else {
    X_global <- NULL
  }
  rm(global_data)
  
  #### Max length of global VAR model ####
  max_lag <- c()
  for (i in countries) {
    max_lag <- append(max_lag, c(model_specs[[i]]$domestic$lags,
                                 model_specs[[i]]$foreign$lags))
  }
  if (global) {
    for (i in countries) {
      max_lag <- append(max_lag, model_specs[[i]]$global$lags)
    }
  }
  max_lag <- max(max_lag)
  for (i in countries){
    model_specs[[i]]$max_lag <- max_lag
  }
  
  #### Generate weight matrices for each country ####
  tv <- !is.na(dim(gvar_weights)[3])
  
  # Check data availability for each country
  n_v <- length(variables)
  var.exists <- as.data.frame(matrix(NA, n_c, n_v))
  row.names(var.exists) <- countries
  names(var.exists) <- variables
  for (i in countries){
    for (j in variables){
      var.exists[i, j] <- is.element(j, index[index[,1] == i, 2])
    }
  }
  k <- dim(index)[1]
  
  W <- c()
  W.names <- c()
  for (i in countries){
    n.endo <- length(index[index[, 1] == i, 1])
    star.vars <- model_specs[[i]]$foreign$variables
    n.star <- length(star.vars)
    if (tv){
      W.i <- array(0, dim = c(n.endo + n.star, k, dim(gvar_weights)[3]))
      dimnames(W.i) <- list(c(index[index[, 1] == i, 2], paste("s.", star.vars, sep = "")),
                            index[, 2], dimnames(gvar_weights)[[3]])
      for (l in 1:dim(gvar_weights)[3]) {
        W.i[1:n.endo, which(index[, 1] == i), l] <- diag(1, n.endo)
        for (j in star.vars){
          w.temp <- gvar_weights[i,,l]
          # Set all values to zero, where variable is not available
          w.temp[!var.exists[, j]] <- 0
          # Extract weights for available values
          exist.weights <- gvar_weights[i, var.exists[, j], l]
          # Make sum over weights of not available values
          non.exist.weights <- sum(gvar_weights[i, !var.exists[, j], l])
          # Recalculate weights
          w.temp[var.exists[, j]] <- exist.weights/(1 - non.exist.weights)
          # Add weights to country's weight matrix
          pos.x <- n.endo + which(star.vars == j)
          pos.y <- which(index[, 2] == j)
          W.i[pos.x, pos.y, l] <- w.temp[var.exists[, j]]
        } 
      }
    } else {
      W.i <- matrix(0, n.endo + n.star, k)
      dimnames(W.i) <- list(c(index[index[, 1] == i, 2], paste("s.", star.vars, sep = "")),
                            index[, 2])
      W.i[1:n.endo, which(index[, 1] == i)] <- diag(1, n.endo)
      for (j in star.vars){
        w.temp <- gvar_weights[i,]
        # Set all values to zero, where variable is not available
        w.temp[!var.exists[, j]] <- 0
        # Extract weights for available values
        exist.weights <- gvar_weights[i, var.exists[, j]]
        # Make sum over weights of not available values
        non.exist.weights <- sum(gvar_weights[i, !var.exists[, j]])
        # Recalculate weights
        w.temp[var.exists[, j]] <- exist.weights/(1 - non.exist.weights)
        # Add weights to country's weight matrix
        pos.x <- n.endo + which(star.vars == j)
        pos.y <- which(index[, 2] == j)
        W.i[pos.x, pos.y] <- w.temp[var.exists[, j]]
      }
    }
    
    W <- c(W, list(W.i))
    W.names <- c(W.names, i)
  }
  names(W) <- W.names
  
  #### Generate vector z = (x, x.star)' for each country ####
  X_temp <- t(X)
  Z <- NULL
  for (j in countries) {
    if (class(W[[j]]) == "matrix"){
      x <- W[[j]] %*% X_temp
    } else {
      x <- matrix(NA, dim(W[[j]])[1], dim(W[[j]])[3])
      for (i in 1:dim(W[[j]])[3]){
        x[, i] <- W[[j]][,, i] %*% X_temp[, i]
      } 
    }
    dimnames(x)[[1]] <- dimnames(W[[j]])[[1]]
    Z <- c(Z, list(x))
    rm(x)
  }
  names(Z) <- countries
  rm(X_temp)
  
  
  #### Produce basic (x, x_star, x_global) objects ####
  X_index <- stats::tsp(X)
  data <- c()
  for (i in 1:length(Z)){
    # Split z into domestic and foreign
    x <- stats::ts(t(matrix(Z[[i]][model_specs[[i]]$domestic$variables,],
                            length(model_specs[[i]]$domestic$variables))),
                   start = X_index[1], frequency = X_index[3], class = c("mts", "ts", "matrix"))
    dimnames(x)[[2]] <- model_specs[[i]]$domestic$variables
    
    x_star <- stats::ts(t(matrix(Z[[i]][-(1:length(model_specs[[i]]$domestic$variables)),],
                                 length(model_specs[[i]]$foreign$variables))),
                        start = X_index[1], frequency = X_index[3], class = c("mts", "ts", "matrix"))
    dimnames(x_star)[[2]] <- model_specs[[i]]$foreign$variables
    
    # Check for proper rank specifications
    if (model_specs[[i]]$type == "VEC" & !any(is.na(model_specs[[i]]$cointegration$rank))) {
      k_temp <- length(model_specs[[i]]$domestic$variables)
      if (length(model_specs[[i]]$cointegration$rank) > (k_temp + 1)) {
        model_specs[[i]]$cointegration$rank <- 0:k_temp
      }
      rm(k_temp)
    }
    
    # Collect data
    temp <- list(data = list(x_domestic = x, x_foreign = x_star),
                 model = model_specs[[i]])
    
    # Check for global variables
    if (!is.null(model_specs[[i]]$global)){
      temp$data$x_global <- X_global
    }
    
    temp$data$weights = W[[i]]
    
    data <- c(data, list(temp))
    rm(list = c("x", "x_star"))
  }
  names(data) <- countries
  
  #### Generate country model for each lag and rank specification ####
  data_temp <- NULL
  names_temp <- NULL
  for (i in countries){
    if (!global) {
      data[[i]]$model$global$lags <- "none"
    }
    if (data[[i]]$model$type == "VAR") {
      data[[i]]$model$cointegration$rank <- "none"
    }
    for (j in data[[i]]$model$domestic$lags) {
      for (k in data[[i]]$model$foreign$lags) {
        for (l in data[[i]]$model$global$lags) {
          for (m in data[[i]]$model$cointegration$rank) {
            temp <- data[[i]]
            temp$model$domestic$lags <- j
            temp$model$foreign$lags <- k
            if (global) {
              temp$model$global$lags <- l
            } else {
              temp$model$global <- NULL
            }
            if (temp$model$type == "VAR") {
              temp$model$cointegration <- NULL
            } else {
              temp$model$cointegration$rank <- m
            }
            
            temp$model$structural <- FALSE
            temp$model$tvp <- FALSE
            temp$model$sv <- FALSE
            temp$model$variable_selection$type <- "none"
            
            data_temp <- c(data_temp, list(temp))
            names_temp <- c(names_temp, i)  
          } 
        }
      }
    }
  }
  
  names(data_temp) <- names_temp
  data <- data_temp
  return(data)
}
