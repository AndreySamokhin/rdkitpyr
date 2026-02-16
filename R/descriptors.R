#==============================================================================#
#' Calculate the exact mass of molecules
#'
#' @description
#'   Calculate the exact mass (monoisotopic mass) for a set of molecules.
#'   Exact mass is computed based on the most abundant isotopes of each element.
#'
#'   The calculation is performed using the \code{ExactMolWt()} function from
#'   the \code{Descriptors} module in RDKit.
#'
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
#' @inheritParams ParseMolecules
#'
#' @return
#'   A numeric vector containing the exact mass for each molecule.
#'   Elements corresponding to invalid molecules are returned as \code{NA}.
#'
#' @examples
#'   # Calculate exact mass for a set of molecules
#'   smiles <- c("CCO", "c1ccccc1", "invalid_molecule")
#'   CalculateExactMass(smiles)
#'   #> 46.04186 78.04695 NA
#'
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
CalculateExactMass <- function(mols,
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


  #--[ Convert to SMILES ]------------------------------------------------------

  py_obj <- reticulate::py$calculate_exact_mass(
    as.list(mols),
    verbose = verbose
  )
  exact_mass <- vapply(reticulate::py_to_r(py_obj), function(a1) {
    if(is.null(a1)) {
      return(NA_real_)
    } else {
      return(a1)
    }
  }, numeric(1L), USE.NAMES = FALSE)
  return(exact_mass)
}



#==============================================================================#
#' Calculate the average molecular weight of molecules
#'
#' @description
#'   Calculate the molecular weight (average mass) for a set of molecules.
#'   Molecular weight is computed based on the natural isotopic abundance of
#'   each element.
#'
#'   The calculation is performed using the \code{MolWt()} function from the
#'   \code{Descriptors} module in RDKit.
#'
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
#' @inheritParams ParseMolecules
#'
#' @return
#'   A numeric vector containing the molecular weight for each molecule.
#'   Elements corresponding to invalid molecules are returned as \code{NA}.
#'
#' @examples
#'   # Calculate exact mass for a set of molecules
#'   smiles <- c("CCO", "c1ccccc1", "invalid_molecule")
#'   CalculateMolecularWeight(smiles)
#'   #> 46.069 78.114 NA
#'
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
CalculateMolecularWeight <- function(mols,
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


  #--[ Convert to SMILES ]------------------------------------------------------

  py_obj <- reticulate::py$calculate_molecular_weight(
    as.list(mols),
    verbose = verbose
  )
  mol_weights <- vapply(reticulate::py_to_r(py_obj), function(a1) {
    if(is.null(a1)) {
      return(NA_real_)
    } else {
      return(a1)
    }
  }, numeric(1L), USE.NAMES = FALSE)
  return(mol_weights)
}


