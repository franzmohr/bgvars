
# Extracts the names of the regressors from a 'ctryvarest' object

.get_regressor_names_var <- function(x) {
  
  names_domestic <- x[["model"]][["domestic"]][["variables"]]
  k_domestic <- length(x[["model"]][["domestic"]][["variables"]])
  p_domestic <- x[["model"]][["domestic"]][["lags"]]
  names_foreign <- x[["model"]][["foreign"]][["variables"]]
  k_foreign <- length(x[["model"]][["foreign"]][["variables"]])
  p_foreign <- x[["model"]][["foreign"]][["lags"]]
  names_global <- x[["model"]][["global"]][["variables"]]
  k_global <- length(x[["model"]][["global"]][["variables"]])
  s_global <- x[["model"]][["global"]][["lags"]]
  if (x[["model"]][["type"]] == "VAR") {
    names_deterministic <- x[["model"]][["deterministic"]]
  }
  if (x[["model"]][["type"]] == "VEC") {
    names_deterministic <- c(x[["model"]][["deterministic"]][["restricted"]],
                             x[["model"]][["deterministic"]][["unrestricted"]])
  }
  x_names <- NULL
  tvp <- x[["model"]][["tvp"]]
  
  if (p_domestic > 0) {
    temp_names <- NULL
    for (i in 1:p_domestic) {
      temp_names <- c(temp_names, paste(names_domestic, ".l", i, sep = ""))
    } 
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(x[["posteriors"]][["foreign"]])) {
    temp_names <- paste0(names_foreign, "*.l", rep(0:p_foreign, each = k_foreign))
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(x[["posteriors"]][["global"]])) {
    temp_names <- paste0(names_global, ".l", rep(0:s_global, each = k_global))
    x_names <- c(x_names, temp_names)
  }
  
  if (!is.null(x[["posteriors"]][["deterministic"]])) {
    x_names <- c(x_names, names_deterministic)
  }
  
  if (!is.null(x[["posteriors"]][["a0"]])) {
    x_names <- c(x_names, names_domestic)
  }
  
  return(x_names)
}

# Extracts the names of the regressors from a 'ctryvecest' object

.get_regressor_names_vec <- function(x) {
  
  names_domestic <- x[["model"]][["domestic"]][["variables"]]
  k_domestic <- length(x[["model"]][["domestic"]][["variables"]])
  p_domestic <- x[["model"]][["domestic"]][["lags"]]
  names_foreign <- x[["model"]][["foreign"]][["variables"]]
  k_foreign <- length(x[["model"]][["foreign"]][["variables"]])
  p_foreign <- x[["model"]][["foreign"]][["lags"]]
  global <- !is.null(x[["model"]][["global"]][["variables"]])
  if (global) {
    names_global <- x[["model"]][["global"]][["variables"]]
    k_global <- length(x[["model"]][["global"]][["variables"]])
    s_global <- x[["model"]][["global"]][["lags"]] 
  }
  names_det_r <- x[["model"]][["deterministic"]][["restricted"]]
  names_det_ur <- x[["model"]][["deterministic"]][["unrestricted"]]
  r <- x[["model"]][["rank"]]
  
  x_names <- NULL
  tvp <- x[["model"]][["tvp"]]
  
  if (p_domestic > 0 & r > 0) {
    x_names <- c(x_names, paste0("l.", names_domestic))
  }
  if (p_foreign > 0 & r > 0) {
    x_names <- c(x_names, paste0("l.", names_foreign, "*"))
  }
  if (global) {
    if (s_global > 0 & r > 0) {
      x_names <- c(x_names, paste0("l.", names_global))
    } 
  }
  if (!is.null(names_det_r) & r > 0) {
    x_names <- c(x_names, names_det_r)
  }
  
  if (!is.null(x[["posteriors"]][["gamma_domestic"]])) {
    if (p_domestic > 1) {
      temp_names <- NULL
      for (i in 1:(p_domestic - 1)) {
        temp_names <- c(temp_names, paste0("d.", names_domestic, ".l", i))
      } 
      x_names <- c(x_names, temp_names)
    }
  }
  
  if (!is.null(x[["posteriors"]][["gamma_foreign"]])) {
    temp_names <- paste0("d.", names_foreign, "*.l", rep(0:(p_foreign - 1), each = k_foreign))
    x_names <- c(x_names, temp_names)
  }
  
  if (global) {
    if (!is.null(x[["posteriors"]][["gamma_global"]])) {
      temp_names <- paste0("d.", names_global, ".l", rep(0:(s_global - 1), each = k_global))
      x_names <- c(x_names, temp_names)
    } 
  }
  
  if (!is.null(x[["posteriors"]][["gamma_deterministic"]])) {
    x_names <- c(x_names, names_det_ur)
  }
  
  if (!is.null(x[["posteriors"]][["a0"]])) {
    x_names <- c(x_names, names_domestic)
  }
  
  return(x_names)
}
