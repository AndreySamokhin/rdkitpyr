SkipIfNoRdkit <- function() {
  if (!.IsRdkitAvailable()) {
    skip("'rdkit' not available for testing.")
  }
}

