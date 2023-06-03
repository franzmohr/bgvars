#' Predict Method for Objects of Class bgvar
#' 
#' Forecasting a Bayesian Global VAR object of class "bgvar" with credible bands.
#' 
#' @param object an object of class "bgvar", usually, a result of a call to \code{\link{combine_models}}.
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
#' @return A time-series object of class "bgvarprd".
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
  
  k <- ncol(object[["data"]][["endogen"]])
  p <- ncol(object[["a"]]) / k^2
  tot <- k * p
  a <- object[["a"]]
  
  model_data <- bvartools::gen_var(object[["data"]][["endogen"]], p = p)
  tt <- nrow(model_data[["data"]][["Y"]])
  
  m <- 0
  if (!is.null(object[["b"]])) {
    a <- cbind(a, object[["b"]])
    m <- ncol(object[["b"]]) / k
    tot <- tot + m
    if (is.null(new_x)) {
      new_x <- matrix(0, n.ahead, m)
    }
    if (NROW(new_x) != n.ahead) {
      stop("Length of argument 'new_x' must be equal to 'n.ahead'.")
    }
  }
  
  if ("c" %in% names(object)) {
    a <- cbind(a, object[["c"]])
    n <- ncol(object[["c"]]) / k
    tot <- tot + n
    if (is.null(new_d)) {
      new_d <- matrix(0, n.ahead, n)
      # Try to find constants and trends automatically
      D_data <- object[["data"]][["deterministic"]]
      const_check <- apply(D_data, 2, function(x){all(x == 1)})
      if (any(const_check)) {
        new_d[, which(const_check)] <- 1
      }
      trend_check <- apply(D_data, 2, function(x){all(x == 1:length(x))})
      if (any(trend_check)) {
        start <- D_data[NROW(D_data), which(trend_check)]
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
  if (m > 0) {
    pos_x <- k * p + 1:m
    pred[pos_x, ] <- 0
    warning("Global variables cannot be used to obtain forecasts yet.")
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
  
  class(result) <- append("bgvarprd", class(result))
  return(result)
}