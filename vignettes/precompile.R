# Pre-compile the expensive vignettes.
#
# The .Rmd files listed below are *generated*: they contain the results of the
# code rather than the code itself, so that `R CMD build` has nothing left to
# evaluate. The real sources are the matching .Rmd.orig files, and this script
# turns the latter into the former. Both the generated .Rmd files and the
# figures under vignettes/figures/ are committed.
#
# Run this after changing a .orig source, or after a change to the package that
# alters vignette output:
#
#   Rscript vignettes/precompile.R          # all vignettes
#   Rscript vignettes/precompile.R GVEC     # just one
#
# Estimating the sub-models takes a while, which is the whole reason this
# script exists. Neither this file nor the .orig sources are shipped in the
# built package -- see .Rbuildignore.

precompile <- function(name) {
  wd <- setwd("vignettes")
  on.exit(setwd(wd), add = TRUE)

  input <- paste0(name, ".Rmd.orig")
  if (!file.exists(input)) {
    stop("no such vignette source: vignettes/", input, call. = FALSE)
  }

  # A chunk that fails must fail the build. knitr's default is to record the
  # error in the output and carry on, which once left a vignette whose every
  # step after a missing folder was an error message.
  knitr::opts_chunk$set(error = FALSE)

  # Every vignette writes to figures/ under its own name, so dropping that
  # prefix leaves no stale plots behind when chunks are renamed or removed, or
  # when one chunk comes to produce a different number of figures than before.
  # This is why the 'fig.path' of a vignette has to start with its own name.
  dir.create("figures", showWarnings = FALSE)
  unlink(list.files("figures", pattern = paste0("^", name, "-"), full.names = TRUE))

  message("Knitting ", input, " ...")
  started <- Sys.time()
  knitr::knit(input, paste0(name, ".Rmd"))
  message(
    "Done in ",
    format(round(difftime(Sys.time(), started), 1)),
    "\n"
  )
}

vignettes <- commandArgs(trailingOnly = TRUE)
if (length(vignettes) == 0) {
  vignettes <- c("GVAR", "GVEC", "TVP-SV-GVAR", "TVP-SV-GVEC",
                 "modelling-submodel")
}

for (vignette in vignettes) {
  precompile(vignette)
}
