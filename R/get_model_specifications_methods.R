
#' @export
#' @method get_model_specifications varxsubmodel
get_model_specifications.varxsubmodel <- function(object, ...) {
  
  result <- data.frame(type = object[["model"]][["type"]],
                       k = object[["model"]][["k_endogen"]],
                       p = object[["model"]][["p_endogen"]],
                       k_exogen = object[["model"]][["k_exogen"]],
                       p_exogen = object[["model"]][["p_exogen"]],
                       m_global = object[["model"]][["m_global"]],
                       s_global = object[["model"]][["s_global"]],
                       n = object[["model"]][["n"]])
  
  # The number of observations is only available, if the object contains data
  tt <- nrow(object[["data"]][["train"]][["y"]])
  if (!is.null(tt)) {
    result[["T"]] <- tt
  }
  
  result[["varsel"]] <- object[["model"]][["varsel"]]

  return(result)
}


#' @export
#' @method get_model_specifications vecxsubmodel
get_model_specifications.vecxsubmodel <- function(object, ...) {
  
  result <- data.frame(type = object[["model"]][["type"]],
                       k = object[["model"]][["k_endogen"]],
                       p = object[["model"]][["p_endogen"]],
                         k_exogen = object[["model"]][["k_exogen"]],
                       p_exogen = object[["model"]][["p_exogen"]],
                       m_global = object[["model"]][["m_global"]],
                       s_global = object[["model"]][["s_global"]],
                       n_unrestricted = object[["model"]][["n"]],
                       n_restricted = object[["model"]][["n_restricted"]],
                       rank = object[["model"]][["rank"]])
  
  # The number of observations is only available, if the object contains data
  tt <- nrow(object[["data"]][["train"]][["y"]])
  if (!is.null(tt)) {
    result[["T"]] <- tt
  }
  
  result[["varsel"]] <- object[["model"]][["varsel"]]
  
  return(result)
}