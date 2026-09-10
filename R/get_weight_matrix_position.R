
.get_weight_matrix_positions <- function(object) {
  
  index <- object[["global"]][["index"]]
  submodels <- unique(index[, "submodel"])
  
  result <- NULL
  for (i in submodels) {
    
    endogen <- object[["submodels"]][[i]][[1]][["model"]][["endogen"]]
    exogen <- object[["submodels"]][[i]][[1]][["model"]][["exogen"]]
    exogen <- gsub(".s", "", exogen)
    
    vars_endogen_old <- index[which(index[, "submodel"] == i) , "variable"]
    n_endogen_old <- length(vars_endogen_old)
    vars_exogen_old <- unique(index[index[, "submodel"] != i, "variable"])
    n_exogen_old <- length(vars_exogen_old)
    n_z_old <- n_endogen_old + n_exogen_old
    
    pos_endogen <- which(index[, "submodel"] == i & index[, "variable"] %in% endogen)
    if (length(pos_endogen) == 0) {
      stop(paste0("For submodel ", submodel, " no variable from argument 'endogen' is available."))
    }
    vars_endogen <- index[pos_endogen, "variable"]
    pos_endogen <- which(vars_endogen_old %in% vars_endogen)
    n_endogen <- length(vars_endogen)
    
    vars_exogen <- unique(index[index[, "submodel"] != i, "variable"])
    pos_exogen <- which(vars_exogen %in% exogen)
    vars_exogen <- vars_exogen[pos_exogen]
    n_exogen <- length(vars_exogen)
    pos_exogen <- n_endogen_old + pos_exogen
    
    n_z <- n_endogen + n_exogen
    
    pos_new <- c(pos_endogen, pos_exogen)
    
    result[[i]] <- list(n_z_old = n_z_old,
                        model = pos_new)
  }
  
  return(result)
}
