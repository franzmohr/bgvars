# Load needed functions from different packages for fast access

#' @importFrom coda thin

# Unload the DLL when the package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("bgvars", libpath)
}