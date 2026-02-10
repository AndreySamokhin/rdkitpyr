#==============================================================================#
#' Infer the type of molecular representation provided as input
#'
#' @param mols
#'   A character vector. Molecular representation, this can be SMILES or InChI
#'   strings.
#'
#' @return
#'   Return a string (\code{inchi}, \code{smiles}, or \code{NA}).
#'
#' @noRd
#==============================================================================#
.InferMolInputType <- function(mols) {
  if (!is.character(mols) || length(mols) == 0L) {
    return(NA_character_)
  }
  if (grepl("^InChI=", mols[[1L]])) {
    return("inchi")
  }
  return("smiles")
}


