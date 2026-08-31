# Unload the DLL when the package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("bgvars", libpath)
}