[![pipeline status](https://git.wur.nl/statistical-genetic-pipeline/statgenQTLxT/badges/master/pipeline.svg)](https://git.wur.nl/statistical-genetic-pipeline/statgenQTLxT/commits/master)
[![coverage report](https://git.wur.nl/statistical-genetic-pipeline/statgenQTLxT/badges/master/coverage.svg)](https://git.wur.nl/statistical-genetic-pipeline/statgenQTLxT/commits/master)

statgenQTLxT
============

R Package for fast multi trait GWAS analysis and genomic prediction

# Implemented functionality

The following functionality has been implemented:

* recoding and imputing markers
* fast multi trait GWAS
* genomic prediction

For single trait GWAS use statgenGWAS (```install.packages("statgenGWAS")```)

# Installation

For direct installation from gitlab use the following code:

``` r
## Replace the location for public and privatekey with your own.
creds <- git2r::cred_ssh_key(publickey = "M:\\.ssh\\id_rsa.pub",
                             privatekey = "M:\\.ssh\\id_rsa")
remotes::install_git(url = "git@git.wur.nl:rossu027/statgenQTLxT.git", credentials = creds)

```
