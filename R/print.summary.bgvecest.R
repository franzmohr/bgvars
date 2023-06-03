#' @include summary.bgvecest.R
#'
#' @export
#' @rdname summary.bgvecest
print.summary.bgvecest <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  for (i in 1:length(x)) {
    if (i != 1) {
      cat("--------------------------------------------------------------------------------\n")
      cat("--------------------------------------------------------------------------------\n\n")
      
    }
    cat("Submodel: ", names(x)[i], "\n")
    print(x[[i]], digits = digits, ...)
  }
}
