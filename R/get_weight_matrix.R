#' Get Weight Matrix
#'
#' Obtains the weight matrix of a sub-model.
#'
#' @param object an object of class 'gvarmodel'.
#'
#' @return A data frame.
#' 
#' @examples
#'
#' @export
get_weight_matrix <- function(object, submodel){
  
  index <- object[["global"]][["index"]]
  vars_endogen <- index[which(index[, "submodel"] == submodel) , "variable"]
  vars_exogen <- unique(index[index[, "submodel"] != submodel, "variable"])
  n_vars <- length(vars_endogen) + length(vars_exogen)
  
  temp <- object[["weights"]][[submodel]][1:n_vars, ]
  dimnames(temp) <- list(c(vars_endogen, vars_exogen),
                         index[, "index"])
  
  return(temp)
}
