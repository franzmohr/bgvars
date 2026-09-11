#' Specifications of the Submodels of a GVAR Model
#'
#' Obtains the model specification of the submodels of a GVAR model.
#'
#' @param object a list of the estimated models of one sub-model, usually one
#' element of \code{object$submodels} after a call to
#' \code{\link[bvartools]{add_posterior_coefficients}}.
#'
#' @return A data frame with one row per model, holding the group the sub-model
#' belongs to, the type of the model, the endogenous, weakly exogenous and global
#' variables with their lag orders, the cointegration rank and the variable
#' selection algorithm. Columns that apply to none of the models -- the rank of a
#' VARX sub-model, for instance -- are dropped.
#'
#' @examples
#'
#' # Load data
#' data("gvar2019")
#' submodel_data <- select_list_elements(gvar2019[["submodel_data"]], c("AT", "US"))
#'
#' # Set up a global model with two candidate lag orders per sub-model
#' object <- create_gvarmodel(submodel_data = submodel_data)
#' object <- add_weight_matrices(object, submodel_data = submodel_data, period = 3)
#' object <- add_submodels(object,
#'                         endogen = c("y", "Dp"), p_endogen = 1:2,
#'                         exogen = c("y", "Dp"), p_exogen = 1,
#'                         iterations = 10, burnin = 10)
#' # Number of iterations and burn-in should be much higher.
#'
#' get_submodel_specifications(object[["submodels"]][["AT"]])
#'
#' @export
get_submodel_specifications <- function(object){
  
  n_models <- length(object)
  
  result <- data.frame(group = rep(NA, n_models),
                       type = rep(NA, n_models),
                       r = rep(NA, n_models),
                       var_domestic = rep(NA, n_models),
                       lag_domestic = rep(NA, n_models),
                       var_foreign = rep(NA, n_models),
                       lag_foreign = rep(NA, n_models),
                       var_global = rep(NA, n_models),
                       lag_global = rep(NA, n_models),
                       varsel = rep(NA, n_models),
                       stringsAsFactors = FALSE)
  
  for (i in 1:n_models) {
    # The models of one sub-model are not named, so there is a group to report
    # only if the caller named them. An all-NA column is dropped below.
    if (!is.null(names(object))) {
      result[i, "group"] <- names(object)[i]
    }
    type <- object[[i]][["model"]][["type"]]
    if (object[[i]][["model"]][["structural"]]) {
     type <- paste0("S", type) 
    }
    if (object[[i]][["model"]][["error"]] %in% c("sv", "sv+covar")) {
      type <- paste0("SV-", type) 
    }
    if (object[[i]][["model"]][["tvp"]]) {
      type <- paste0("TVP-", type) 
    }
    result[i, "type"] <- type
    rm(type)
    result[i, "var_domestic"] <- paste(object[[i]][["model"]][["endogen"]], collapse = ", ")
    result[i, "lag_domestic"] <- object[[i]][["model"]][["p_endogen"]]
    result[i, "var_foreign"] <- paste(object[[i]][["model"]][["exogen"]], collapse = ", ")
    result[i, "lag_foreign"] <- object[[i]][["model"]][["p_exogen"]]
    if (object[[i]][["model"]][["m_global"]] > 0) {
      result[i, "var_global"] <- paste(object[[i]][["model"]][["global"]], collapse = ", ")
      result[i, "lag_global"] <- object[[i]][["model"]][["s_global"]]
    }
    if (!is.null(object[[i]][["model"]][["rank"]])) {
      result[i, "r"] <- object[[i]][["model"]][["rank"]] 
    }
    if (!is.null(object[[i]][["model"]][["varsel"]])) {
      result[i, "varsel"] <- object[[i]][["model"]][["varsel"]] 
    }
  }
  
  result <- result[, !unlist(lapply(result, function(x) {all(is.na(x))}))]
  
  return(result)
}
