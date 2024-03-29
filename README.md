
<!-- README.md is generated from README.Rmd. Please edit that file -->

# statgenQTLxT

[![](https://www.r-pkg.org/badges/version/statgenQTLxT)](https://www.r-pkg.org/pkg/statgenQTLxT)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/statgenQTLxT)](https://www.r-pkg.org/pkg/statgenQTLxT)
[![R-CMD-check](https://github.com/Biometris/statgenQTLxT/workflows/R-CMD-check/badge.svg)](https://github.com/Biometris/statgenQTLxT/actions?workflow=R-CMD-check)
[![codecov](https://codecov.io/gh/Biometris/statgenQTLxT/branch/main/graph/badge.svg)](https://app.codecov.io/gh/Biometris/statgenQTLxT)

**statgenQTLxT** is an R package for fast multi-trait and multi-trial
Genome Wide Association Studies (GWAS).

statgenQTLxT has extensive options for summarizing and visualizing
results. It builds on the `statgenGWAS` package (for single trait GWAS)
which is available from
[CRAN](https://biometris.github.io/statgenGWAS/). The package uses data
structures and plots defined in the `statgenGWAS` package. It is
recommended to read the vignette of this package, accessible in R via
`vignette(package = "statgenGWAS")` or online at
<https://biometris.github.io/statgenGWAS/articles/GWAS.html> to get a
general idea of those.

For a full description of the theoretical background and a fully worked
example see the
[**vignette**](https://biometris.github.io/statgenQTLxT/articles/statgenQTLxT.html).

## Installation

- Install from CRAN:

``` r
install.packages("statgenQTLxT")
```

- Install latest development version from GitHub (requires
  [remotes](https://github.com/r-lib/remotes) package):

``` r
remotes::install_github("Biometris/statgenQTLxT", ref = "develop", dependencies = TRUE)
```
