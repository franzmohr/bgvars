#' Specifications of the Submodels of a GVAR Model
#'
#' Obtains the model specification of the submodels of a GVAR model.
#'
#' @param object an object of class 'submodelestlist', usually,
#' a result of a call to \code{\link[bvartools]{draw_posterior}}.
#'
#' @return A data frame.
#' 
#' @export
# The example this function used to carry was written for the API that
# preceded the rename to 'endogen'/'exogen': it called create_specifications(),
# create_models() and draw_posterior(), none of which exist any more. The
# function itself has not been migrated either -- see
# tests/testthat/test-known-issues.R -- so there is no example that would run.
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
    result[i, "group"] <- names(object)[i]
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
    result[i, "var_domestic"] <- paste(object[[i]][["model"]][["domestic_vars"]], collapse = ", ")
    result[i, "lag_domestic"] <- object[[i]][["model"]][["p_domestic"]]
    result[i, "var_foreign"] <- paste(object[[i]][["model"]][["foreign_vars"]], collapse = ", ")
    result[i, "lag_foreign"] <- object[[i]][["model"]][["p_foreign"]]
    if (object[[i]][["model"]][["m"]] > 0) {
      result[i, "var_global"] <- paste(object[[i]][["model"]][["global_vars"]], collapse = ", ")
      result[i, "lag_global"] <- object[[i]][["model"]][["s"]] 
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
