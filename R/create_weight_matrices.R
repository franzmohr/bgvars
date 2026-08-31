
.create_weight_matrices <- function(submodel, index, weights) {
  
  weights <- weights[[submodel]]
  vars_weights <- dimnames(weights)[[2]]
  
  pos_endogen <- which(index[, "submodel"] == submodel)
  vars_endogen <- index[pos_endogen, "variable"]
  n_endogen <- length(vars_endogen)
  
  vars_exogen <- unique(index[index[, "submodel"] != submodel, "variable"])
  n_exogen <- length(vars_exogen)
  
  n_z <- n_endogen + n_exogen
  n_model <- nrow(index)
  tt <- nrow(weights)
  
  w <- matrix(0, n_z * tt, n_model)
  
  for (i in 1:tt) {
    
    # Diagonal matris corresponding to a sub-model's own endogenous variables
    w[(i - 1) * n_z + 1:n_endogen, pos_endogen] <- diag(1, n_endogen)
    
    # Weights for weakly endogenous variables
    for (var_i in vars_exogen) {
      index_i <- index[index[, "variable"] == var_i & index[, "submodel"] != submodel, ]
      pos_x <- which(vars_exogen == var_i)
      pos_exogen <- index_i[, "id"]
      pos_weights <- which(vars_weights %in% index_i[, "submodel"])
      w[(i - 1) * n_z + n_endogen + pos_x, pos_exogen] <- weights[i, pos_weights] / sum(weights[i, pos_weights])
    }
  }
  
  return(w)
}
