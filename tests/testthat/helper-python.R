SkipIfNoRdkit <- function() {
  if (!.IsRdkitAvailable(initialize = TRUE)) {
    skip("'rdkit' not available for testing.")
  }
}

