
<!-- README.md is generated from README.Rmd. Please edit that file -->

# statgenQTLxT

[![](https://www.r-pkg.org/badges/version/statgenQTLxT)](https://www.r-pkg.org/pkg/statgenQTLxT)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/statgenQTLxT)](https://www.r-pkg.org/pkg/statgenQTLxT)

**statgenQTLxT** is an R package for fast multi trait Genome Wide
Association Studies (GWAS).

statgenQTLxT has extensive options for summarizing and visualizing
results. It builds on the `statgenGWAS` package (for single trait GWAS)
which is available from
[CRAN](https://biometris.github.io/statgenGWAS/). The package uses data
structures and plots defined in the `statgenGWAS` package. It is
recommended to read the vignette of this package, accessible in R via
`vignette(package = "statgenGWAS")` or online at
<https://biometris.github.io/statgenGWAS/articles/GWAS.html> to get a
general idea of those.

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
