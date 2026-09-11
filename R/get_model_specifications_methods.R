
#' Get Model Specifications
#'
#' Obtains the specifications of a sub-model of a global model.
#'
#' @param object an object of class \code{'varxsubmodel'}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A one-row data frame with the specifications of the sub-model: the
#' number of endogenous variables and their lag order, the number of weakly
#' exogenous and of global variables and their lag orders, the number of
#' deterministic terms, the number of observations if the object
#' contains data, and the variable selection algorithm.
#'
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


#' Get Model Specifications
#'
#' Obtains the specifications of a sub-model of a global model.
#'
#' @param object an object of class \code{'vecxsubmodel'}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A one-row data frame with the specifications of the sub-model: the
#' number of endogenous variables and their lag order, the number of weakly
#' exogenous and of global variables and their lag orders, the number of
#' deterministic terms, restricted and unrestricted, and the rank of the
#' cointegration matrix, the number of observations if the object
#' contains data, and the variable selection algorithm.
#'
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