[![pipeline status](https://git.wur.nl/rossu027/genStatPipeline/badges/master/pipeline.svg)](https://git.wur.nl/rossu027/genStatPipeline/commits/master)

# genStatPipeline

R Genetic Statistical Pipeline for 
* recoding and imputing markers
* single trait GWAS
* multi trait GWAS

# Installation

For direct installation from gitlab use the following code:

``` r
## Replace the location for public and privatekey with your own.
creds <- git2r::cred_ssh_key(publickey = "M:\\.ssh\\id_rsa.pub",
                             privatekey = "M:\\.ssh\\id_rsa")
devtools::install_git(url = "git@git.wur.nl:rossu027/statGenPipeline.git",
                      credentials = creds)

```
