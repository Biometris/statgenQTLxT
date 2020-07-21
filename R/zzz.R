.onLoad <- function(libname = find.package("statgenQTLxT"),
                    pkgname = "statgenQTLxT"){
  ## CRAN Note avoidance.
  if (getRversion() >= "2.15.1") {
    utils::globalVariables("j")
  }
  invisible()
}

.onUnload <- function(libpath) {
  library.dynam.unload("statgenQTLxT", libpath)
}
