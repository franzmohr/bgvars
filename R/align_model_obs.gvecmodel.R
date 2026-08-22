#' Align Observations across Models
#' 
#' Restricts each sub-model in an object of class 'gvecmodel' to the set of
#' observations common to all sub-models, ensuring that comparisons are computed
#' on the same underlying sample.
#' 
#' @param object an object of class 'gvecmodel'.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return An object of class 'gvecmodel'.
#' 
#' @export
#' @method align_model_obs gvecmodel
align_model_obs.gvecmodel <- function(object, ...) {
  
  # Get sample sizes
  avail <- NULL
  for (i in 1:length(object[["submodels"]])) {
    if ("modellist" %in% class(object[["submodels"]][[i]])) {
      for (j in 1:length(object[["submodels"]][[i]])) {
        result_i <- stats::tsp(object[["submodels"]][[i]][[j]][["data"]][["train"]][["y"]])
        avail <- rbind(avail, result_i)
      }
    } else {
      result_i <- stats::tsp(object[["submodels"]][[i]][["data"]][["train"]][["y"]])
      avail <- rbind(avail, result_i)
    }
  }

  min_date <- max(avail[, 1]) # First period overall is latest in start
  max_date <- min(avail[, 2]) # First period overall is earlierst in end
  
  object <- stats::window(x = object, start = min_date, end = max_date)
  
}

