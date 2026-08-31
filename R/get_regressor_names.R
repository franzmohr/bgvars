
# Extracts the names of the regressors from a 'varxsubmodelest' object

.get_regressor_names_var <- function(x) {
  
  names_domestic <- x[["model"]][["domestic_vars"]]
  k_domestic <- x[["model"]][["k_domestic"]]
  p_domestic <- x[["model"]][["p_domestic"]]
  names_foreign <- x[["model"]][["foreign_vars"]]
  k_foreign <- x[["model"]][["k_foreign"]]
  p_foreign <- x[["model"]][["p_foreign"]]
  names_global <- x[["model"]][["global_vars"]]
  m <- x[["model"]][["m"]]
  s <- x[["model"]][["s"]]
  global <- m > 0
  n <- x[["model"]][["n"]]
  names_deterministic <- x[["model"]][["deterministic"]]
  
  x_names <- NULL
  tvp <- x[["model"]][["tvp"]]
  
  if (p_domestic > 0) {
    temp_names <- NULL
    for (i in 1:p_domestic) {
      temp_names <- c(temp_names, paste(names_domestic, ".l", i, sep = ""))
    } 
    x_names <- c(x_names, temp_names)
  }
  
  temp_names <- paste0(names_foreign, "*.l", rep(0:p_foreign, each = k_foreign))
  x_names <- c(x_names, temp_names)
  
  if (m > 0) {
    temp_names <- paste0(names_global, ".l", rep(0:s, each = k_global))
    x_names <- c(x_names, temp_names)
  }
  
  if (n > 0) {
    x_names <- c(x_names, names_deterministic)
  }
  
  if (x[["model"]][["structural"]]) {
    x_names <- c(x_names, names_domestic)
  }
  
  return(x_names)
}

# Extracts the names of the regressors from a 'vecxsubmodelest' object

.get_regressor_names_vec <- function(x) {
  
  names_domestic <- x[["model"]][["domestic_vars"]]
  k_domestic <- x[["model"]][["k_domestic"]]
  p_domestic <- x[["model"]][["p_domestic"]]
  names_foreign <- x[["model"]][["foreign_vars"]]
  k_foreign <- x[["model"]][["k_foreign"]]
  p_foreign <- x[["model"]][["p_foreign"]]
  names_global <- x[["model"]][["global_vars"]]
  m <- x[["model"]][["m"]]
  s <- x[["model"]][["s"]]
  global <- m > 0
  r <- x[["model"]][["rank"]]
  n_restricted <- x[["model"]][["n_restricted"]]
  names_det_r <- dimnames(x[["data"]][["determ_restricted"]])[[2]]
  n_unrestricted <- x[["model"]][["n_unrestricted"]]
  names_det_ur <- dimnames(x[["data"]][["determ_unrestricted"]])[[2]]
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
    if (s > 0 & r > 0) {
      x_names <- c(x_names, paste0("l.", names_global))
    } 
  }
  if (n_unrestricted > 0 & r > 0) {
    x_names <- c(x_names, names_det_r)
  }
  
  if (p_domestic > 1) {
    temp_names <- NULL
    for (i in 1:(p_domestic - 1)) {
      temp_names <- c(temp_names, paste0("d.", names_domestic, ".l", i))
    } 
    x_names <- c(x_names, temp_names)
  }
  
  temp_names <- paste0("d.", names_foreign, "*.l", rep(0:(p_foreign - 1), each = k_foreign))
  x_names <- c(x_names, temp_names)
  
  
  if (global) {
    temp_names <- paste0("d.", names_global, ".l", rep(0:(s_global - 1), each = k_global))
    x_names <- c(x_names, temp_names)
  }
  
  if (n_unrestricted > 0) {
    x_names <- c(x_names, names_det_ur)
  }
  
  if (x[["model"]][["structural"]]) {
    x_names <- c(x_names, names_domestic)
  }
  
  return(x_names)
}
