.onLoad <- function(libname = find.package("genStatPipeline"),
                    pkgname = "genStatPipeline"){
  ## CRAN Note avoidance.
  if (getRversion() >= "2.15.1") {
    utils::globalVariables("j")
  }
  invisible()
}

.onUnload <- function(libpath) {
  library.dynam.unload("genStatPipeline", libpath)
}
