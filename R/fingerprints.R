#==============================================================================#
#' Calculate MACCS fingerprints for a set of molecules
#'
#' @description
#'   Calculate MACCS (Molecular ACCess System) fingerprints for a set of
#'   molecules. Each fingerprint is a fixed-length binary vector representing
#'   the presence or absence of predefined structural features.
#'
#'   The fingerprints are calculated using the \code{GetMACCSKeysFingerprint()}
#'   function from the \code{rdMolDescriptors} module in RDKit.
#'
#'   Invalid molecules are represented by \code{NA} vectors. Row order is
#'   preserved so that the output aligns with the input.
#'
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
#' @inheritParams ParseMolecules
#'
#' @return
#'   A matrix of integers (0/1) with one row per molecule and one column per
#'   MACCS key (length 167). Rows corresponding to invalid molecules contain
#'   \code{NA}.
#'
#'   Additionally, the \code{"valid"} attribute is attached to indicate which
#'   molecules were successfully processed.
#'
#' @examples
#'   # Calculate MACCS fingerprints
#'   smiles <- c("CCO", "c1ccccc1", "invalid_molecule")
#'   fps <- CalculateMaccsFingerprints(smiles)
#'
#'   # Get the number of fingerprints (columns)
#'   ncol(fps)
#'   #> 167
#'
#'   # Check which molecules were successfully processed
#'   attr(fps, "valid")
#'   #> TRUE TRUE FALSE
#'
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
CalculateMaccsFingerprints <- function(mols,
                                       verbose = FALSE) {

  .EnsurePythonReady()


  #--[ Check input arguments ]--------------------------------------------------

  # 'mols'
  if ((!is.character(mols) || length(mols) == 0L) &&
      (!is.list(mols) || length(mols) == 0L ||
       (!is.null(mols[[1L]]) &&
        !inherits(mols[[1L]], "rdkit.Chem.rdchem.Mol")))) {
    stop("'mols' must be a character vector (SMILES or InChI) ",
         "or a list of RDKit Mol objects.")
  }

  # 'verbose'
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a logical value.")
  }


  #--[ Calculate fingerprints ]-------------------------------------------------

  # rdkit.Chem.MACCSkeys.GenMACCSKeys()
  # https://rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html
  py_obj <- reticulate::py$calculate_maccs_fps(
    as.list(mols),
    verbose = verbose
  )
  fp_size <- 167L
  fps <- t(vapply(reticulate::py_to_r(py_obj), function(a1) {
    if(is.null(a1)) {
      return(rep(NA_integer_, fp_size))
    } else {
      return(a1)
    }
  }, integer(fp_size), USE.NAMES = FALSE))
  attr(fps, "valid") <- !is.na(fps[, 1L])
  return(fps)
}


