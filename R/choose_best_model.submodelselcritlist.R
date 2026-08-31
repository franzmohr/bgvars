#' Choose Best Model
#' 
#' Chooses the best model according the selection criteria in an object of class
#' 'submodelselcritlist'.
#' 
#' @param object object of class 'submodelselcritlist', usually, a result of a call
#' to \code{\link[bvartools]{selection_criteria}}.
#' @param criterion the selection criterion that should be plotted. Available choices
#' are \code{"LL"}, \code{"AIC"}, \code{"BIC"} (default), \code{"HQ"}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details
#' If argument \code{criterion} is \code{"LL"}, the model with the maximum value is chosen,
#' otherwise, the model with the minimum value.
#' 
#' @returns A vector of integers containing the positions of the best models
#' in the list provided in argument \code{object}.
#' 
#' @export
choose_best_model.submodelselcritlist <- function(object, criterion = "BIC", ...) {
  
  res <- lapply(object, function(y) {y[["summary"]]})
  res <- do.call("rbind", res)
  rownames(res) <- NULL
  res <- cbind(data.frame(submodel = names(object)), res)
  
  submodels <- unique(names(object))
  
  pos_final_models <- c()
  for (i in submodels) {
    pos_overall <- which(names(object) == i)
    if (criterion == "LL") {
      pos_best_model_i <- which(res[res$submodel == i, criterion] == max(res[res$submodel == i, criterion]))
    } else {
      pos_best_model_i <- which(res[res$submodel == i, criterion] == min(res[res$submodel == i, criterion])) 
    }
    pos_final_i <- pos_overall[pos_best_model_i]
    pos_final_models <- append(pos_final_models, pos_final_i)
  }
  
  return(pos_final_models)
}
