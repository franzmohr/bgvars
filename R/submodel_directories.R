#' The Directories of an Export
#'
#' The directory each sub-model of an exported global model lies in, which is
#' what an external program that works on the files is pointed at.
#'
#' @param object an object of class 'gvarexport', from
#' \code{\link{read_export}}.
#' @param submodels character vector of the sub-models to return the
#' directories of, or \code{NULL}, the default, for all of them.
#'
#' @details
#' \code{\link[bayests]{write_to_hdf5}} gives every sub-model a directory of
#' its own and every model of it a file in that directory, so a sub-model's
#' directory is a natural unit of work: it holds the models of one country and
#' nothing else.
#'
#' What this is for is handing those directories to something that does the work
#' on the files rather than in the session. The posterior is the largest thing a
#' model has, and \code{\link{map_submodels}}, which reads a model, gives it to a
#' function and writes it back, is the wrong shape for producing one: it would
#' copy the draws into R and out again for nothing. With
#' \code{\link[bayests]{bayests_files}} a global model is estimated as
#' \preformatted{run <- bayests_files(executable = "/opt/bayests/bin/bayests")
#' run(submodel_directories(export), jobs = 6)}
#' which starts six BayesTS processes, each working through one country's models
#' where they lie. The draws are the same as those of
#' \code{\link[bayests]{add_posterior_coefficients}}, since a model is drawn
#' with the seed in its file.
#'
#' Note that the parallelism is in the processes rather than in a cluster of R
#' workers. An idle R worker is not cheap -- on Windows a fresh one commits about
#' 2 GB before it does anything and about 4 GB once it has loaded a package of
#' this kind -- so a cluster that only waits for external processes can cost more
#' memory than the processes themselves.
#'
#' @return A named character vector of directories, one per sub-model, in the
#' order of the manifest.
#'
#' @examples
#'
#' data("gvar2023")
#' submodel_data <- select_list_elements(gvar2023[["submodel_data"]],
#'                                       c("AT", "DE", "US"))
#' object <- create_gvarmodel(submodel_data = submodel_data,
#'                            global_data = gvar2023[["global_data"]])
#' object <- add_weight_matrices(object = object,
#'                               submodel_data = submodel_data, period = 3)
#' object <- add_submodels(object, endogen = c("y", "Dp"), p_endogen = 1,
#'                         exogen = c("y", "Dp"), p_exogen = 1,
#'                         deterministic = "const",
#'                         iterations = 20, burnin = 10)
#'
#' folder <- file.path(tempdir(), "bgvars-example-directories")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
#'
#' export <- read_export(folder)
#' submodel_directories(export)
#' submodel_directories(export, submodels = "US")
#'
#' @family export
#' @export
submodel_directories <- function(object, submodels = NULL) {

  if (!inherits(object, "gvarexport")) {
    stop("Argument 'object' must be of class 'gvarexport'. Use read_export().")
  }

  manifest <- object[["manifest"]]
  names_of <- unique(manifest[, "submodel"])
  if (!is.null(submodels)) {
    unknown <- setdiff(submodels, names_of)
    if (length(unknown) > 0) {
      stop("This export does not contain the sub-model(s) ",
           paste0(unknown, collapse = ", "), ".")
    }
    names_of <- names_of[names_of %in% submodels]
  }

  vapply(names_of, function(submodel) {
    rows <- manifest[, "submodel"] == submodel
    directory <- unique(dirname(manifest[rows, "file"]))
    if (length(directory) != 1) {
      stop("The models of sub-model ", submodel, " lie in ", length(directory),
           " directories, so there is no one directory to name.")
    }
    file.path(object[["folder"]], directory)
  }, character(1))
}
