#' Predict Method for Objects of Class bgvar
#' 
#' Forecasting a Bayesian Global VAR object of class \code{"bgvar"} with credible bands.
#' 
#' @param object an object of class \code{"bgvar"}, usually, a result of a call to \code{\link{combine_submodels}}.
#' @param n.ahead number of steps ahead at which to predict.
#' @param new_x a matrix of new non-deterministic, exogenous variables. Must have \code{n.ahead} rows.
#' @param new_d a matrix of new deterministic variables. Must have \code{n.ahead} rows.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param ... additional arguments.
#' 
#' @details The function produces \code{n.ahead} forecasts for the GVAR model
#' \deqn{y_t = \sum_{l = 1}^{p} G_{l} y_{t-i} + \sum_{l = 0}^{s} H_{l} x_{t-i} + D d_t + G^{-1}_{0} u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' @return A time-series object of class \code{"bgvarprd"}.
#' 
#' @examples
#' # Load data
#' data("gvar2016")
#' 
#' # Create regions
#' temp <- create_regions(country_data = gvar2016$country_data,
#'              weight_data = gvar2016$weight_data,
#'              region_weights = gvar2016$region_weights,
#'              regions = list(EA =  c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")),
#'              period = 3)
#' 
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' global_data = gvar2016$global_data
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
#' # Overwrite country-specific specifications
#' model_specs[["US"]][["domestic"]][["variables"]] <- c("y", "Dp", "r")
#' model_specs[["US"]][["foreign"]][["variables"]] <- c("y", "Dp")
#' 
#' # Create estimation objects
#' country_models <- create_models(country_data = country_data,
#'                                 weight_data = weight_data,
#'                                 global_data = global_data,
#'                                 model_specs = model_specs)
#' 
#' # Add priors
#' models_with_priors <- add_priors(country_models,
#'                                  coef = list(v_i = 1 / 9, v_i_det = 1 / 10),
#'                                  sigma = list(df = 3, scale = .0001))
#' 
#' # Obtain posterior draws
#' object <- draw_posterior(models_with_priors)
#' 
#' # Solve GVAR
#' gvar <- combine_submodels(object)
#' 
#' # Obtain forecasts
#' gvar_prd <- predict(gvar)
#' 
#' # Plot forecast
#' plot(gvar_prd, variable = c("US", "y")) 
#' 
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
#' 
#' 
#' @export
predict.bgvar <- function(object, ..., n.ahead = 10, new_x = NULL, new_d = NULL, ci = .95) {
  
  # Dev specs
  # n.ahead = 10; new_x = NULL; new_d = rep(1, 10); ci = .95
  
  k <- length(object[["model"]][["endogen"]][["variables"]]) # Endogenous variables
  p <- object[["model"]][["endogen"]][["lags"]]  # Determine number of lags
  tot <- k * p # Total number of endogenous regressors
  
  # Check if global
  global <- !is.null(object[["model"]][["global"]])
  s <- NULL
  if (global) {
    k_glo <- length(object[["model"]][["global"]][["variables"]]) # Endogenous variables
    s <- object[["model"]][["global"]][["lags"]]  # Determine number of lags
  }
  
  # Generate a simple VAR model
  model_data <- bvartools::gen_var(data = object[["data"]][["endogen"]], p = p,
                                   exogen = object[["data"]][["global"]], s = s)
  tt <- nrow(model_data[["data"]][["Y"]]) # Number of observations
  
  # Endogenous variables
  a <- object[["a"]]
  
  # Global variables
  m <- 0
  if (global) {
    a <- cbind(a, object[["b"]])
    n_glo <- k_glo * (s + 1)
    tot <- tot + n_glo
    if (is.null(new_x)) {
      warning("Argument 'new_x' empty. Replacing it by a matrix of zeros.")
      new_x <- matrix(0, n.ahead, n_glo)
    }
    if (NROW(new_x) != n.ahead) {
      stop("Length of argument 'new_x' must be equal to 'n.ahead'.")
    }
    new_x <- rbind(model_data[["data"]][["Z"]][tt, k * p + 1:n_glo], new_x)
  }
  
  # Determinisitc terms
  if ("c" %in% names(object)) {
    a <- cbind(a, object[["c"]])
    n <- ncol(object[["c"]]) / k
    tot <- tot + n
    if (is.null(new_d)) {
      new_d <- matrix(0, n.ahead, n)
      # Try to find constants and trends automatically
      D_data <- object[["data"]][["deterministic"]]
      # Check if any column of deterministic data contains only ones -> const
      const_check <- apply(D_data, 2, function(x){all(x == 1)})
      if (any(const_check)) {
        # Fill respective column with ones
        new_d[, which(const_check)] <- 1
      }
      # Check if any column of deterministic data contains a linear sequence -> trenc
      trend_check <- apply(D_data, 2, function(x){all(x == 1:length(x))})
      if (any(trend_check)) {
        # Determine last value of trend
        start <- D_data[NROW(D_data), which(trend_check)]
        # Add trend to respective column with correct starting period
        new_d[, which(trend_check)] <- start:(start + n.ahead - 1)
      }
    }
    if (NROW(new_d) != n.ahead) {
      stop("Length of argument 'new_d' must be equal to 'n.ahead'.")
    }
  }
  
  pos_y <- 1:(k * p)
  # Generate empty prediction matrix
  pred <- matrix(NA, tot, n.ahead + 1)
  
  # Add starting values
  pred[pos_y, 1] <- t(model_data[["data"]][["Y"]][tt:(tt - p + 1),])
  if (global) {
    pos_x <- k * p + 1:n_glo
    pred[pos_x, ] <- t(new_x)
  }
  if ("c" %in% names(object)) {
    pos_d <- (tot - n + 1):tot
    pred[pos_d, ] <- cbind(object[["data"]][["deterministic"]][tt, ], t(new_d))
  }
  
  # Calculate forecasts
  draws <- nrow(a)
  result <- array(NA, dim = c(k, n.ahead, draws))
  pb <- utils::txtProgressBar(style = 3)
  for (draw in 1:draws) {
    for (i in 1:n.ahead) {
      
      a0_i <- solve(matrix(object[["a0"]][draw, ], k))
      
      # Generate random error
      temp <- eigen(matrix(object[["sigma"]][draw, ], k))
      u <- temp$vectors %*% diag(sqrt(temp$values), k) %*% t(temp$vectors) %*% stats::rnorm(k)
      
      pred[1:k, i + 1] <- a0_i %*% matrix(a[draw, ], k) %*% pred[, i] + a0_i %*% u
      if (p > 1) {
        for (j in 1:(p - 1)) {
          pred[j * k + 1:k, i + 1] <- pred[(j - 1) * k + 1:k, i]
        } 
      }
    }
    result[,, draw] <- pred[1:k, -1]
    utils::setTxtProgressBar(pb, draw / draws)
  }
  close(pb)
  
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  temp <- apply(result, c(2, 1) , stats::quantile, probs = c(ci_low, .5, ci_high)) 
  result <- c()
  for (i in 1:k) {
    result <- c(result, list(stats::ts(t(temp[,, i]),
                                       start = stats::tsp(object[["data"]][["endogen"]])[2],
                                       frequency = stats::tsp(object[["data"]][["endogen"]])[3])))
    names(result)[i] <- paste(object[["index"]][i, ], collapse = "_")
  }
  
  result <- list("y" = model_data[["data"]][["Y"]],
                 "fcst" = result,
                 "index" = object[["index"]])
  
  class(result) <- c("bgvarprd", "list")
  return(result)
}