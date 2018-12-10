## Create data for testing functions.
## Restricted and anonymized version of DROPs data.
load(system.file("extdata", "testData.RData", package = "genStatPipeline"))

# Export to package
devtools::use_data(Y, K, X, overwrite = TRUE)
