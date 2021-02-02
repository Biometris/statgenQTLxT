pkgsUpdate <- function(quiet = TRUE,
                       instPkgdown = FALSE) {
  rver <- paste0(R.version$major, ".", R.version$minor)
  instld <- as.data.frame(installed.packages())
  upd <- rownames(instld[instld$Built != rver, ])
  install.packages(upd, quiet = FALSE)
  remotes::install_deps(dependencies = TRUE, quiet = quiet)
  cat("INSTALLED:\n")
  instld <- as.data.frame(installed.packages())
  rownames(instld) <- NULL
  print(instld[, c("Package", "Version")])
  return(invisible(TRUE))
}
