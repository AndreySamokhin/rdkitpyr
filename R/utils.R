#==============================================================================#
#' Infer the type of molecular representation provided as input
#'
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
#'
#' @return
#'   Return a string (\code{"inchi"}, \code{"smiles"}, \code{"rdkit_mol"} or
#'   \code{NA}).
#'
#' @noRd
#==============================================================================#
.InferMolInputType <- function(mols) {
  if ((!is.character(mols) || length(mols) == 0L) &&
      (!is.list(mols) || !inherits(mols[[1L]], "rdkit.Chem.rdchem.Mol"))) {
    return(NA_character_)
  }
  if (grepl("^InChI=", mols[[1L]])) {
    return("inchi")
  } else if (inherits(mols[[1L]], "rdkit.Chem.rdchem.Mol")) {
    return("rdkit_mol")
  }
  return("smiles")
}


