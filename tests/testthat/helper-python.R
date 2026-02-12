SkipIfNoRdkit <- function() {
  if (!reticulate::py_module_available("rdkit")) {
    skip("'rdkit' not available for testing.")
  }
}

