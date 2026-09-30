#' Apply a Function to the Models of an Export
#'
#' Reads every model of an exported global model, applies a function to it,
#' writes it back and drops it, so that a step over a global model costs the
#' memory of one model per worker rather than of the whole model.
#'
#' @param object an object of class 'gvarexport', from
#' \code{\link{open_export}}.
#' @param f a function taking an estimable model -- a 'bvarmodel' or a
#' 'bvecmodel' -- and returning one, or, with \code{write = FALSE}, returning
#' anything at all. If it has an argument named \code{index}, the position of the
#' model in the manifest is passed to it, which is what numbers seeds and what
#' lets a caller tell the models apart.
#' @param ... further arguments passed to \code{f}.
#' @param write whether what \code{f} returns is written back to the model's
#' file. \code{TRUE}, the default, is a step of an estimation, and \code{f} must
#' then return a model. \code{FALSE} is a step that reads: \code{f} may return
#' anything, the export is left exactly as it is, and the results come back in a
#' named list instead of the object. It is what a diagnostic over an export
#' uses, as \code{write} of \code{\link[bvartools]{map_models}} does for a folder
#' of single models.
#' @param submodels character vector of the sub-models the step is applied to,
#' or \code{NULL}, the default, for all of them. The models of the other
#' sub-models are left untouched, and so are their rows of the manifest.
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
#' With \code{write = FALSE} nothing is written and the manifest is not rebuilt,
#' so a diagnostic cannot change what the export is. The results are named by
#' the sub-model and the file that holds the model, as in \code{AT/001}, because
#' a sub-model has one model per candidate specification and the file is what
#' tells them apart. \code{\link[bvartools]{time_variation_test}} for an export
#' is such a step.
#'
#' A step over a subset is a step over the whole export applied to part of it,
#' not a step over a smaller export: the position a model has in the whole
#' manifest is what is passed as \code{index} and what numbers its seed, so a
#' run taken sub-model by sub-model draws what a run over all of them draws.
#' That is what makes a long estimation resumable -- the sub-models that are
#' already done are passed over rather than read and written again -- and it is
#' how a specification too large to estimate in one sitting is estimated in
#' several. The estimation steps for an export pass \code{submodels} on, so
#' \code{add_posterior_coefficients(export, submodels = "AT")} simulates the
#' models of one country.
#'
#' With more than one core, \code{f} is sent to the workers as it stands, and
#' the workers load the installed packages. A function that calls something a
#' package does not export, or something only a newer version of it has, fails
#' there rather than in the session it was written in. Pass what such a function
#' needs as an argument instead.
#'
#' @return With \code{write = TRUE}, the object in \code{object}, with its
#' manifest brought up to date, invisibly. With \code{write = FALSE}, a named
#' list of what \code{f} returned for each model.
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
#' # One sub-model at a time, which is how a long estimation is resumed
#' export <- map_submodels(export, add_initial_values, submodels = "AT")
#'
#' # A step that reads rather than writes
#' map_submodels(export, function(model) nrow(model$data$train$y),
#'               write = FALSE)
#'
#' @family export
#' @export
map_submodels <- function(object, f, ..., submodels = NULL, write = TRUE,
                          cores = 1, export = NULL) {

  if (!inherits(object, "gvarexport")) {
    stop("Argument 'object' must be of class 'gvarexport'. Use open_export().")
  }
  if (!is.logical(write) || length(write) != 1 || is.na(write)) {
    stop("Argument 'write' must be TRUE or FALSE.")
  }
  f <- match.fun(f)

  manifest <- object[["manifest"]]
  # The rows the step is applied to. Their positions in the whole manifest are
  # kept, since that is what a model's seed is numbered by.
  rows <- seq_len(nrow(manifest))
  if (!is.null(submodels)) {
    unknown <- setdiff(submodels, manifest[, "submodel"])
    if (length(unknown) > 0) {
      stop("This export does not contain the sub-model(s) ",
           paste0(unknown, collapse = ", "), ".")
    }
    rows <- which(manifest[, "submodel"] %in% submodels)
  }
  if (length(rows) == 0) {
    return(if (write) invisible(object) else list())
  }

  # Everything a model needs travels in its own task -- the function and the
  # arguments included -- so that nothing is passed alongside the tasks and the
  # worker function takes one argument. Sending arguments through
  # parallel::parLapplyLB() next to them makes the call depend on how that
  # function forwards its dots.
  arguments <- list(...)
  tasks <- lapply(rows, function(i) {
    list(file = file.path(object[["folder"]], manifest[i, "file"]), index = i,
         f = f, arguments = arguments, write = write)
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

  if (!write) {
    names(specs) <- paste(manifest[rows, "submodel"],
                          sub("\\.h5$", "", basename(manifest[rows, "file"])),
                          sep = "/")
    return(specs)
  }

  for (j in seq_along(specs)) {
    i <- rows[j]
    manifest[i, ] <- .manifest_row(specs[[j]], manifest[i, "submodel"],
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

  # A step that reads hands back whatever the function made of the model, and
  # the model itself is dropped with the worker.
  if (!isTRUE(task[["write"]])) {
    return(model)
  }

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
