# Solving a global model from its export, a piece of the chain at a time.
#
# Stacking the sub-models is the one step that needs every sub-model at once,
# but it needs them one draw at a time: a draw of the global model is built from
# that draw of every sub-model and from nothing else. So the draws are read in
# pieces, each piece is stacked with the code that stacks a model in memory, and
# the pieces are put end to end. What a solve then holds of the sub-models is
# one piece rather than the whole chain, which for a global model of time
# varying sub-models is the difference between a few gigabytes and tens of them.
#
# The global model's own draws are held whole. They are the smaller half -- one
# k by k block per lag and draw, against a coefficient path per period and draw
# in every time varying sub-model -- and every method that works on the solved
# model expects them together.

.solve_export <- function(object, period = NULL, submodels = NULL, chunk = 100) {

  if (!inherits(object, "gvarexport")) {
    stop("Argument 'object' must be of class 'gvarexport'. Use open_export().")
  }
  chunk <- max(1L, as.integer(chunk))

  manifest <- object[["manifest"]]
  if (!is.null(submodels)) {
    manifest <- .select_manifest_rows(manifest, submodels)
  }
  counts <- table(manifest[, "submodel"])
  if (any(counts > 1)) {
    stop("The export holds more than one model for ",
         paste(names(counts)[counts > 1], collapse = ", "),
         ". Name the model to take per sub-model in 'submodels', as ",
         "read_gvar_from_folder() does.")
  }

  # Before a file is read: the manifest records the algorithm, and a closed form
  # has no draws to stack. .combine_submodels() would say the same on the first
  # chunk, after reading it.
  algorithm <- .discount_submodel_algorithm(manifest[, "algorithm"])
  if (!is.null(algorithm)) {
    .refuse_discount_submodel(algorithm)
  }

  chain <- .export_chain(object[["folder"]], manifest)
  pieces <- split(seq_len(chain[["draws"]]),
                  ceiling(seq_len(chain[["draws"]]) / chunk))

  message("Solving the global model from ", object[["folder"]], ", ",
          chain[["draws"]], " draws in ", length(pieces),
          if (length(pieces) == 1) " piece ..." else " pieces ...")

  solved <- NULL
  for (piece in pieces) {
    part <- read_gvar_from_folder(object[["folder"]], submodels = submodels,
                                  draws = piece)
    solved <- .stack_solved(solved, suppressMessages(.combine_submodels(part, period = period)))
    rm(part)
  }

  # The draws of the pieces were read as chains of their own, so the labels of
  # the chain the export holds are put back here, and the number of draws the
  # specification reports is the whole chain rather than the last piece.
  as_draws <- function(x) {
    coda::mcmc(x, start = chain[["start"]], end = chain[["end"]],
               thin = chain[["thin"]])
  }
  solved[["posterior"]][["a"]][["coeffs"]] <- as_draws(solved[["posterior"]][["a"]][["coeffs"]])
  solved[["posterior"]][["u_sigma_inv"]][["coeffs"]] <-
    as_draws(solved[["posterior"]][["u_sigma_inv"]][["coeffs"]])
  solved[["g"]] <- as_draws(solved[["g"]])
  solved[["model"]][["iterations"]] <- as.integer(chain[["draws"]])

  solved
}

# The chain of an export: how many draws the sub-models hold and how they are
# labelled. Read from the first file, since a global model whose sub-models
# disagree on either cannot be stacked and .combine_submodels() says so.
.export_chain <- function(folder, manifest) {

  filename <- file.path(folder, manifest[1, "file"])
  file <- hdf5r::h5file(filename, mode = "r")
  on.exit(if (file$is_valid) file$close_all(), add = TRUE)

  if (!file$exists("posterior")) {
    stop("Sub-model ", manifest[1, "submodel"], " has not been estimated yet, ",
         "so the global model cannot be solved.")
  }
  dataset <- file[["posterior"]][["u_sigma_inv"]][["coeffs"]]
  draws <- dataset$dims[1]
  labels <- c(start = hdf5r::h5attr(dataset, "start"),
              end = hdf5r::h5attr(dataset, "end"),
              thin = hdf5r::h5attr(dataset, "thin"))
  file$close_all()

  list(draws = as.integer(draws), start = labels[["start"]],
       end = labels[["end"]], thin = labels[["thin"]])
}

# Two solved models, one after the other. Everything but the draws describes the
# model rather than the piece it was built from, so it is taken from the first.
.stack_solved <- function(first, second) {

  if (is.null(first)) {
    return(second)
  }

  rows <- function(name) {
    rbind(unclass(first[["posterior"]][[name]][["coeffs"]]),
          unclass(second[["posterior"]][[name]][["coeffs"]]))
  }
  first[["posterior"]][["a"]][["coeffs"]] <- rows("a")
  first[["posterior"]][["u_sigma_inv"]][["coeffs"]] <- rows("u_sigma_inv")
  first[["g"]] <- rbind(unclass(first[["g"]]), unclass(second[["g"]]))

  first
}
