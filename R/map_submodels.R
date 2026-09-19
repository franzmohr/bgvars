#' Apply a Function to the Models of an Export
#'
#' Reads every model of an exported global model, applies a function to it,
#' writes it back and drops it, so that a step over a global model costs the
#' memory of one model per worker rather than of the whole model.
#'
#' @param object an object of class 'gvarexport', from
#' \code{\link{open_export}}.
#' @param f a function taking an estimable model -- a 'bvarmodel' or a
#' 'bvecmodel' -- and returning one. If it has an argument named \code{index},
#' the position of the model in the manifest is passed to it, which is what
#' numbers seeds and what lets a caller tell the models apart.
#' @param ... further arguments passed to \code{f}.
#' @param cores the number of worker processes the models are handled on.
#' Defaults to one. With more than one the models are spread over a socket
#' cluster; every model has its own file, so the workers share nothing.
#' @param export character vector of names of objects in the calling
#' environment that the workers need, for instance the sampler that
#' \code{posterior_function} refers to. Ignored with one core.
#'
#' @details
#' The unit is one model, which is what a file of the export holds. A model is
#' read with \code{\link[bvartools]{read_model_from_hdf5}}, passed to \code{f},
#' and written back to its own file through a temporary file that replaces the
#' old one once it is complete, so an interrupted step leaves either the model
#' as it was or the model as \code{f} made it, never half of it.
#'
#' The manifest is rebuilt from the specifications of the models that come back
#' and written to \code{model.h5}, so a step that changes what a model is -- its
#' rank, its lag orders, the chain it asks for -- leaves the table describing the
#' files. It records what a model is, not how many draws are left in it, so
#' thinning the draws does not change it.
#'
#' The methods of the estimation steps for an export call this, so a pipeline
#' rarely calls it directly. It is what to use for a step the package has no
#' method for.
#'
#' With more than one core, \code{f} is sent to the workers as it stands, and
#' the workers load the installed packages. A function that calls something a
#' package does not export, or something only a newer version of it has, fails
#' there rather than in the session it was written in. Pass what such a function
#' needs as an argument instead.
#'
#' @return The object in \code{object}, with its manifest brought up to date,
#' invisibly.
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
#' folder <- file.path(tempdir(), "bgvars-example-map")
#' dir.create(folder, showWarnings = FALSE)
#' write_to_hdf5(object, folder = folder, overwrite = TRUE, mc.cores = 1)
#'
#' export <- open_export(folder)
#'
#' # What add_priors() for an export does
#' export <- map_submodels(export, function(model) {
#'   add_priors(model, coef = list(v_i = 1), sigma = list(df = 3, scale = 1e-8))
#' })
#'
#' @family export
#' @export
map_submodels <- function(object, f, ..., cores = 1, export = NULL) {

  if (!inherits(object, "gvarexport")) {
    stop("Argument 'object' must be of class 'gvarexport'. Use open_export().")
  }
  f <- match.fun(f)

  manifest <- object[["manifest"]]
  # Everything a model needs travels in its own task -- the function and the
  # arguments included -- so that nothing is passed alongside the tasks and the
  # worker function takes one argument. Sending arguments through
  # parallel::parLapplyLB() next to them makes the call depend on how that
  # function forwards its dots.
  arguments <- list(...)
  tasks <- lapply(seq_len(nrow(manifest)), function(i) {
    list(file = file.path(object[["folder"]], manifest[i, "file"]), index = i,
         f = f, arguments = arguments)
  })

  cores <- max(1L, as.integer(cores))
  cores <- min(cores, length(tasks))

  specs <- if (cores < 2L) {
    lapply(tasks, .map_one_model)
  } else {
    cluster <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterEvalQ(cluster, suppressMessages(library(bgvars)))
    if (!is.null(export)) {
      parallel::clusterExport(cluster, export, envir = parent.frame())
    }
    parallel::parLapplyLB(cluster, tasks, .map_one_model)
  }

  for (i in seq_along(specs)) {
    manifest[i, ] <- .manifest_row(specs[[i]], manifest[i, "submodel"],
                                   manifest[i, "file"])
  }
  object[["manifest"]] <- manifest
  .write_export_manifest(object)

  invisible(object)
}

# One model: read, apply, write back, and return its specification for the
# manifest. The model itself never leaves the worker.
.map_one_model <- function(task) {

  f <- task[["f"]]
  model <- bvartools::read_model_from_hdf5(filename = task[["file"]])

  arguments <- c(list(model), task[["arguments"]])
  if ("index" %in% names(formals(f))) {
    arguments <- c(arguments, list(index = task[["index"]]))
  }
  model <- do.call(f, arguments)

  if (!inherits(model, "bvarmodel") && !inherits(model, "bvecmodel")) {
    stop("The function returned an object of class ",
         paste(class(model), collapse = ", "), " for ", task[["file"]],
         " rather than a model.")
  }

  # Through a temporary file, so that an interrupted write cannot leave a model
  # that is neither the old one nor the new one. bvartools' writer refuses an
  # existing file, which is why the target is removed rather than written over.
  temporary <- paste0(task[["file"]], ".writing")
  if (file.exists(temporary)) {
    unlink(temporary)
  }
  bvartools::write_to_hdf5(model, filename = temporary)
  unlink(task[["file"]])
  if (!file.rename(temporary, task[["file"]])) {
    stop("The model of ", task[["file"]], " was written to ", temporary,
         " and could not be moved into place.")
  }

  model[["model"]]
}
