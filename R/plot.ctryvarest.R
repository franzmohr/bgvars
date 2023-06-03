#' Plotting Draws of a VARX Submodel of a GVAR Model
#' 
#' A plot function for objects of class \code{"ctryvarest"} for visual inspection
#' of posterior draws.
#' 
#' @param x an object of class \code{"ctryvarest"}, usually, a result of a call to \code{\link{draw_posterior}}.
#' @param ci interval used to calculate credible bands for time-varying parameters.
#' @param type either \code{"hist"} (default) for histograms, \code{"trace"} for a trace plot,
#' or \code{"boxplot"} for a boxplot. Only used for parameter draws of constant coefficients.
#' @param variables character vector of variables that should be plotted. Default is \code{"all"}.
#' Other options are \code{"domestic"}, \code{"foreign"}, \code{"global"}, \code{"deterministic"}
#' and \code{"sigma"}.
#' @param ctry character (optional). Name of the country, which will be shown in the title of the output.
#' @param ... further graphical parameters.
#' 
#' @export
plot.ctryvarest <- function(x, ci = 0.95, type = "hist", variables = "all", ctry = NULL, ...) {
  
  if (!type %in% c("hist", "trace", "boxplot")) {
    stop("Argument 'type' must be 'hist', 'trace' or 'boxplot'.")
  }
  
  if (!variables %in% c("all", "domestic", "foreign", "deterministic", "sigma")) {
    stop("Invalid specification of argument 'variables'.")
  }
  
  names_domestic <- x[["model"]][["domestic"]][["variables"]]
  k_domestic <- length(x[["model"]][["domestic"]][["variables"]])
  p_domestic <- x[["model"]][["domestic"]][["lags"]]
  names_foreign <- x[["model"]][["foreign"]][["variables"]]
  k_foreign <- length(x[["model"]][["foreign"]][["variables"]])
  p_foreign <- x[["model"]][["foreign"]][["lags"]]
  names_global <- x[["model"]][["global"]][["variables"]]
  global <- !is.null(x[["model"]][["global"]][["variables"]])
  if (global) {
    names_global <- x[["model"]][["global"]][["variables"]]
    k_global <- length(x[["model"]][["global"]][["variables"]])
    s_global <- x[["model"]][["global"]][["lags"]] 
  } else {
    k_global <- 0
    s_global <- 0
  }
  names_deterministic <- x[["model"]][["deterministic"]]
  if (is.null(names_deterministic)) {
    k_deterministic <- 0
  } else {
    k_deterministic <- length(x[["model"]][["deterministic"]]) 
  }
  
  tt <- nrow(x[["data"]][["Y"]])
  tsp_info <- stats::tsp(x[["data"]][["Y"]])
  structural <- x[["model"]][["structural"]]
  ci_low <- (1 - ci) / 2
  ci_high <- 1 - ci_low
  y_names <- dimnames(x[["data"]][["domestic"]])[[2]]
  x_names <- .get_regressor_names_var(x)
  lab_size <- .05
  mar_orig <- graphics::par("mar")
  
  # Calculate number of columns in final layout
  n_tot <- 0
  if ("domestic" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_domestic * p_domestic
  }
  if ("foreign" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_foreign * (p_foreign + 1)
  }
  if ("global" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_global * (s_global + 1)
  }
  if ("deterministic" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_deterministic
  }
  if (!is.null(x[["posteriors"]][["a0"]])) {
    n_tot <- n_tot + k_domestic
  }
  if ("sigma" %in% names(x[["posteriors"]])) {
    n_tot <- n_tot + k_domestic
  }
  
  if (variables == "domestic") {
    if (p_domestic > 0) {
      n_tot <- k_domestic * p_domestic
      x_names <- x_names[1:n_tot]
    }
  }
  if (variables == "foreign") {
    n_tot <- k_foreign * (p_foreign + 1)
    x_names <- x_names[k_domestic * p_domestic + 1:n_tot]
  }
  if (variables == "global") {
    if (global) {
      n_tot <- k_global * (s_global + 1)
      x_names <- x_names[k_domestic * p_domestic + k_foreign * p_foreign + 1:n_tot] 
    }
  }
  if (variables == "deterministic") {
    if (k_deterministic > 0) {
      n_tot <- k_deterministic
      x_names <- x_names[k_domestic * p_domestic + k_foreign * p_foreign + k_global * s_global + 1:n_tot] 
    }
  }
  if (variables == "sigma") {
    n_tot <- k_domestic
  }
  
  # Define layout
  mat <- matrix(NA_integer_,
                k_domestic + 2 , # k_domestic + title row + column names
                n_tot + 1) # n_tot + row names
  mat[1, ] <- 1
  mat[-1, 1] <- c(0, 2:(k_domestic + 1))
  mat[2, -1] <- (k_domestic + 1) + 1:n_tot
  mat[-(1:2), -1] <- matrix(1:(k_domestic * n_tot) + k_domestic + n_tot + 1, k_domestic, n_tot)
  graphics::layout(mat,
                   widths = c(lab_size, rep((1 - lab_size) / n_tot, n_tot)),
                   heights = c(.07, lab_size, rep((1 - lab_size) / k_domestic, k_domestic)))
  
  # Title
  title_text <- "Bayesian "
  tvp <- x[["model"]][["tvp"]]
  if (tvp) {
    title_text <- paste0(title_text, "TVP-")
  }
  if (x[["model"]][["sv"]]) {
    title_text <- paste0(title_text, "SV-")
  }
  if (x[["model"]][["structural"]]) {
    title_text <- paste0(title_text, "S")
  }
  title_text <- paste0(title_text, "VARX submodel")
  if (!is.null(ctry)) {
    title_text <- paste0(title_text, " - ", ctry)
  }
  
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new(); graphics::text(0.5, 0.5, labels = title_text, cex = 1.5)
  # Fill rows
  graphics::par(mar = c(3, 0, 0, 0))
  for (j in y_names) {
    graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
  }
  # Fill columns
  graphics::par(mar = c(0, 0, 0, 0))
  if (variables %in% c("all", "domestic", "foreign", "global", "deterministic")) {
    for (j in x_names) {
      graphics::plot.new(); graphics::text(0.5, 0.5, labels = j, adj = 0.5)
    }
  }
  if (variables %in% c("all", "sigma")) {
    for (j in y_names) {
      graphics::plot.new(); graphics::text(0.5, 0.5, labels = paste0("Sigma\n", j), adj = 0.5)
    } 
  }
  
  graphics::par(mar = c(3, 2.1, .5, 1))
  
  if (variables %in% c("all", "domestic")) {
    if (!is.null(x[["posteriors"]][["domestic"]])) {
      for (i in 1:p_domestic) {
        var_pos <- ((i - 1) * k_domestic * k_domestic) + 1:(k_domestic * k_domestic)
        if (tvp) {
          for (j in var_pos) {
            draws <- .tvpribbon(x[["A"]], j, ci_low, ci_high)
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "")
          }
        } else {
          for (j in var_pos) {
            if (all(x[["posteriors"]][["domestic"]][, j] == x[["posteriors"]][["domestic"]][1, j])) {
              graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["domestic"]][1, j], adj = 0.5)
            } else {
              if (type == "hist") {
                graphics::hist(x[["posteriors"]][["domestic"]][, j], plot = TRUE, main = NA)  
              }
              if (type == "trace") {
                stats::ts.plot(x[["posteriors"]][["domestic"]][, j], xlab = "")
              }
              if (type == "boxplot") {
                graphics::boxplot(x[["posteriors"]][["domestic"]][, j])
              }
            } 
          }
        }
      }
    }
  }
  
  if (variables %in% c("all", "foreign")) {
    if (!is.null(x[["posteriors"]][["foreign"]])) {
      for (i in 1:(p_foreign + 1)) {
        var_pos <- ((i - 1) * k_domestic * k_foreign) + 1:(k_domestic * k_foreign)
        if (tvp) {
          for (j in var_pos) {
            draws <- .tvpribbon(x[["B"]], j, ci_low, ci_high)
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "")
          }
        } else {
          for (j in var_pos) {
            if (all(x[["posteriors"]][["foreign"]][, j] == x[["posteriors"]][["foreign"]][1, j])) {
              graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["foreign"]][1, j], adj = 0.5)
            } else {
              if (type == "hist") {
                graphics::hist(x[["posteriors"]][["foreign"]][, j], plot = TRUE, main = NA)  
              }
              if (type == "trace") {
                stats::ts.plot(x[["posteriors"]][["foreign"]][, j], xlab = "")
              }
              if (type == "boxplot") {
                graphics::boxplot(x[["posteriors"]][["foreign"]][, j])
              }
            }
          }
        }
      }
    }
  }
  
  if (variables %in% c("all", "global")) {
    if (!is.null(x[["posteriors"]][["global"]])) {
      for (i in 1:(s_global + 1)) {
        var_pos <- ((i - 1) * k_domestic * k_global) + 1:(k_domestic * k_global)
        if (tvp) {
          for (j in var_pos) {
            draws <- .tvpribbon(x[["B"]], j, ci_low, ci_high)
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "")
          }
        } else {
          for (j in var_pos) {
            if (all(x[["posteriors"]][["global"]][, j] == x[["posteriors"]][["global"]][1, j])) {
              graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["global"]][1, j], adj = 0.5)
            } else {
              if (type == "hist") {
                graphics::hist(x[["posteriors"]][["global"]][, j], plot = TRUE, main = NA)  
              }
              if (type == "trace") {
                stats::ts.plot(x[["posteriors"]][["global"]][, j], xlab = "")
              }
              if (type == "boxplot") {
                graphics::boxplot(x[["posteriors"]][["global"]][, j])
              }
            }
          }
        }
      }
    }
  }
  
  if (variables %in% c("all", "deterministic")) {
    if (!is.null(x[["posteriors"]][["deterministic"]])) {
      if (tvp) {
        for (j in 1:NCOL(x[["posteriors"]][["deterministic"]][[1]])) {
          draws <- .tvpribbon(x[["posteriors"]][["deterministic"]], j, ci_low, ci_high)
          stats::tsp(draws) <- tsp_info
          stats::ts.plot(draws, xlab = "")
        }
      } else {
        for (j in 1:NCOL(x[["posteriors"]][["deterministic"]])) {
          if (all(x[["posteriors"]][["deterministic"]][, j] == x[["posteriors"]][["deterministic"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["deterministic"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["deterministic"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["deterministic"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["deterministic"]][, j])
            }
          }
        }
      }
    }
  }
  
  if (variables %in% c("all")) {
    if (!is.null(x[["posteriors"]][["a0"]])) {
      if (tvp) {
        for (j in 1:NCOL(x[["posteriors"]][["a0"]][[1]])) {
          draws <- .tvpribbon(x[["posteriors"]][["a0"]], j, ci_low, ci_high)
          if (all(draws[, 1] == draws[1, 1])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
          } else {
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "")
          }
        }
      } else {
        for (j in 1:NCOL(x[["posteriors"]][["a0"]])) {
          if (all(x[["posteriors"]][["a0"]][, j] == x[["posteriors"]][["a0"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["a0"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["a0"]][, j], plot = TRUE, main = NA)   
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["a0"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["a0"]][, j])
            }
          }
        }
      }
    }
  }
  
  if (variables %in% c("all", "sigma")) {
    if (!is.null(x[["posteriors"]][["sigma"]])) {
      var_pos <- 1:(k_domestic * k_domestic)
      if (x[["model"]][["sv"]]) {
        for (j in var_pos) {
          draws <- .tvpribbon(x[["posteriors"]][["sigma"]], j, ci_low, ci_high)
          if (all(draws[, 1] == draws[1, 1])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = draws[1, 2], adj = 0.5)
          } else {
            stats::tsp(draws) <- tsp_info
            stats::ts.plot(draws, xlab = "") 
          }
        }
      } else {
        for (j in var_pos) {
          if (all(x[["posteriors"]][["sigma"]][, j] == x[["posteriors"]][["sigma"]][1, j])) {
            graphics::plot.new(); graphics::text(0.5, 0.5, labels = x[["posteriors"]][["sigma"]][1, j], adj = 0.5)
          } else {
            if (type == "hist") {
              graphics::hist(x[["posteriors"]][["sigma"]][, j], plot = TRUE, main = NA)  
            }
            if (type == "trace") {
              stats::ts.plot(x[["posteriors"]][["sigma"]][, j], xlab = "")
            }
            if (type == "boxplot") {
              graphics::boxplot(x[["posteriors"]][["sigma"]][, j])
            }
          }
        }
      }
    }
  }
  
  graphics::par(mar = mar_orig)
  graphics::layout(matrix(1))
}


