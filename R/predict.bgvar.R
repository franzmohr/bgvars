#' Predict Method for Objects of Class bgvar
#' 
#' Forecasting a Bayesian Global VAR object of class "bgvar" with credible bands.
#' 
#' @param object an object of class "bgvar", usually, a result of a call to \code{\link{solve_gvar}}.
#' @param n.ahead number of steps ahead at which to predict.
#' @param new_x a matrix of new non-deterministic, exogenous variables. Must have \code{n.ahead} rows.
#' @param new_D a matrix of new deterministic variables. Must have \code{n.ahead} rows.
#' @param ci a numeric between 0 and 1 specifying the probability mass covered by the
#' credible intervals. Defaults to 0.95.
#' @param ... additional arguments.
#' 
#' @details The function produces \code{n.ahead} forecasts for the GVAR model
#' \deqn{y_t = \sum_{l = 1}^{p} G_{l} y_{t-i} + \sum_{l = 0}^{s} H_{l} x_{t-i} + D d_t + G^{-1}_{0} u_t,}
#' with \eqn{u_t \sim N(0, \Sigma)}.
#' 
#' @return A time-series object of class "bgvarprd".
#' 
#' @references
#' 
#' Lütkepohl, H. (2007). \emph{New introduction to multiple time series analyis}. Berlin: Springer.
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
#' # Generate EA area region with 2 year rolling window weights
#' ea <- c("AT", "BE", "DE", "ES", "FI", "FR", "IT", "NL")
#' temp <- create_regions(country_data = country_data,
#'                        regions = list("EA" = ea),
#'                        period = 2,
#'                        region_weights = region_weights,
#'                        weight_data = weight_data)
#' country_data <- temp$country_data
#' weight_data <- temp$weight_data
#' 
#' # Generate weight matrices as 2 year, rolling window averages
#' gvar_weights <- create_weights(weight_data = weight_data, period = 2,
#'                                country_data = country_data)
#' 
#' # Create an object with country model specifications
#' model_specs <- create_specifications(country_data = country_data,
#'                                      global_data = global_data,
#'                                      variables = c("y", "Dp", "r"),
#'                                      countries = c("EA", "US", "JP", "CA", "MX", "GB"),
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
#' gvar_est <- estimate_gvar(object, iterations = 100, burnin = 10, thin = 2)
#' # Note that the number of iterations and burnin should be much higher.
#' 
#' # Solve GVAR
#' gvar_solved <- solve_gvar(gvar_est)
#' 
#' # Forecasts
#' gvar_predict <- predict(gvar_solved)
#'
#' # Plot Forecasts
#' plot(gvar_predict, variable = c("US", "Dp"))
#' 
#' @export
predict.bgvar <- function(object, ..., n.ahead = 10, new_x = NULL, new_D = NULL, ci = .95) {
  k <- ncol(object$data$X)
  p <- ncol(object$g) / k^2
  tot <- k * p
  g <- object$g
  
  model_data <- bvartools::gen_var(object$data$X, p = p)
  t <- ncol(model_data$Y)
  
  if (!is.null(object$h)) {
    g <- cbind(g, object$h)
    m <- ncol(object$h) / k
    tot <- tot + m
    if (is.null(new_x)) {
      new_x <- matrix(0, n.ahead, m)
    }
    if (NROW(new_x) != n.ahead) {
      stop("Length of argument 'new_x' must be equal to 'n.ahead'.")
    }
  }
  
  if ("d" %in% names(object)) {
    g <- cbind(g, object$d)
    n <- ncol(object$d) / k
    tot <- tot + n
    if (is.null(new_D)) {
      new_D <- matrix(0, n.ahead, n)
      # Try to find constants and trends automatically
      D_data <- object$data$deterministics  
      const_check <- apply(D_data, 2, function(x){all(x == 1)})
      if (any(const_check)) {
        new_D[, which(const_check)] <- 1
      }
      trend_check <- apply(D_data, 2, function(x){all(x == 1:length(x))})
      if (any(trend_check)) {
        start <- D_data[NROW(D_data), which(trend_check)]
        new_D[, which(trend_check)] <- start:(start + n.ahead - 1)
      }
    }
    if (NROW(new_D) != n.ahead) {
      stop("Length of argument 'new_D' must be equal to 'n.ahead'.")
    }
  }
  
  pos_y <- 1:(k * p)
  pred <- matrix(NA, tot, n.ahead + 1)
  pred[pos_y, 1] <- model_data$Y[, t:(t - p + 1)]
  if (!is.null(new_x)) {
    pos_x <- k * p + 1:m
    pred[pos_x, ] <- cbind(object$data$exogen[t, ], t(new_x))
  }
  if ("d" %in% names(object)) {
    pos_d <- (tot - n + 1):tot
    pred[pos_d, ] <- cbind(object$data$deterministics[t, ], t(new_D))
  }
  
  draws <- nrow(g)
  result <- array(NA, dim = c(k, n.ahead, draws))
  pb <- utils::txtProgressBar(style = 3)
  for (draw in 1:draws) {
    for (i in 1:n.ahead) {
      g0_i <- matrix(object$g0_i[draw, ], k)
      temp <- eigen(matrix(object$sigma[draw, ], k))
      u <- temp$vectors %*% diag(sqrt(temp$values), k) %*% t(temp$vectors) %*% stats::rnorm(k)
      pred[1:k, i + 1] <- matrix(g[draw, ], k) %*% pred[, i] + g0_i %*% u
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
    result <- c(result, list(stats::ts(t(temp[,, i]), start = stats::tsp(object$data$X)[2],
                                       frequency = stats::tsp(object$data$X)[3])))
    names(result)[i] <- paste(object$index[i, ], collapse = "_")
  }
  
  result <- list("y" = model_data$Y,
                 "fcst" = result,
                 "index" = object$index)
  
  class(result) <- append("bgvarprd", class(result))
  return(result)
}