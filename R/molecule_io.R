#==============================================================================#
#' Convert molecular representations to SMILES strings
#'
#' @description
#'   Convert molecular representations to SMILES strings.
#'
#'   The input must be a character vector of either SMILES or InChI identifiers;
#'   mixing the two in the same vector is not allowed.
#'
#'   SMILES inputs can also be used to transform their representation, for
#'   example, to canonical SMILES or to include stereochemistry.
#'
#' @param mols
#'   A character vector. Molecular representation, this can be SMILES or InChI
#'   strings. All elements must be of the same type.
#' @param isomeric
#'   A logical value. If \code{TRUE}, include stereochemistry in the output.
#' @param kekule
#'   A logical value. If \code{TRUE}, kekulize aromatic bonds.
#' @param explicit_bonds
#'   A logical value. If \code{TRUE}, make all bonds explicit.
#' @param explicit_hydrogens
#'   A logical value. If \code{TRUE}, make all hydrogens explicit.
#' @param canonical
#'   A logical value. If \code{TRUE}, generate canonical SMILES.
#'
#' @return
#'   A character vector. SMILES strings. Elements that cannot be converted are
#'   returned as \code{NA}.
#'
#' @examples
#'   # Convert a vector of InChI identifiers to canonical SMILES
#'   inchi <- c("InChI=1S/C2H6/c1-2/h1-2H3",
#'              "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")
#'   ConvertToSmiles(inchi)
#'   #> "CC"
#'   #> "c1ccccc1"
#'
#'   # Convert a vector of SMILES to SMILES with Kekulized aromatic bonds
#'   smiles <- c("c1ccccc1", "c1ccc2ccccc2c1")
#'   ConvertToSmiles(smiles, kekule = TRUE)
#'   #> "C1=CC=CC=C1"
#'   #> "C1=CC=C2C=CC=CC2=C1"
#'
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
ConvertToSmiles <- function(mols,
                            isomeric = TRUE,
                            kekule = FALSE,
                            explicit_bonds = FALSE,
                            explicit_hydrogens = FALSE,
                            canonical = TRUE) {

  #--[ Pre-processing ]---------------------------------------------------------

  .EnsurePythonReady()
  mol_input_type <- .InferMolInputType(mols[1L])


  #--[ Check input arguments ]--------------------------------------------------

  # 'mols'
  if (is.na(mol_input_type)) {
    stop("'mols' must be a vector of SMILES or InChI identifiers.")
  }

  # 'isomeric'
  if (!is.logical(isomeric) || length(isomeric) != 1L) {
    stop("'isomeric' must be a logical value.")
  }

  # 'kekule'
  if (!is.logical(kekule) || length(kekule) != 1L) {
    stop("'kekule' must be a logical value.")
  }

  # 'explicit_bonds'
  if (!is.logical(explicit_bonds) || length(explicit_bonds) != 1L) {
    stop("'explicit_bonds' must be a logical value.")
  }

  # 'explicit_hydrogens'
  if (!is.logical(explicit_hydrogens) || length(explicit_hydrogens) != 1L) {
    stop("'explicit_hydrogens' must be a logical value.")
  }

  # 'canonical'
  if (!is.logical(canonical) || length(canonical) != 1L) {
    stop("'canonical' must be a logical value.")
  }


  #--[ Convert to SMILES ]------------------------------------------------------

  # A single element R vector is converted to a Python scalar. To overcome this
  # behavior R vectors are represented as lists explicitly.

  py_obj <- reticulate::py$convert_to_smiles(
    as.list(mols),
    isomericSmiles = isomeric,
    kekuleSmiles = kekule,
    allBondsExplicit = explicit_bonds,
    allHsExplicit = explicit_hydrogens,
    canonical = canonical
  )
  smiles <- vapply(reticulate::py_to_r(py_obj), function(a1) {
    if(is.null(a1)) {
      return(NA_character_)
    } else {
      return(a1)
    }
  }, character(1L), USE.NAMES = FALSE)
  return(smiles)
}



#==============================================================================#
#' Convert molecular representations to InChI strings
#'
#' @description
#'   Convert molecular representations to InChI strings.
#'
#'   InChI identifiers can technically be provided as input. In this case, the
#'   output is expected to be identical to the input. This can be useful to test
#'   RDKit's consistency with challenging molecules.
#'
#' @param mols
#'   A character vector. Molecular representation, this can be SMILES or InChI
#'   strings. All elements must be of the same type.
#'
#' @return
#'   A character vector. InChI strings. Elements that cannot be converted are
#'   returned as \code{NA}.
#'
#' @examples
#'   # Convert a vector of SMILES to InChI identifiers
#'   smiles <- c("CC", "CCC")
#'   ConvertToInchi(smiles)
#'   #> "InChI=1S/C2H6/c1-2/h1-2H3"
#'   #> "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3"
#'
#'   # Providing InChI as input returns identical output
#'   ConvertToInchi("InChI=1S/CH4/h1H4")
#'   #> "InChI=1S/CH4/h1H4"
#'
#'
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
ConvertToInchi <- function(mols) {

  #--[ Pre-processing ]---------------------------------------------------------

  .EnsurePythonReady()
  mol_input_type <- .InferMolInputType(mols[1L])


  #--[ Check input arguments ]--------------------------------------------------

  if (is.na(mol_input_type)) {
    stop("'mols' must be a vector of SMILES or InChI identifiers.")
  }


  #--[ Convert Smiles ]---------------------------------------------------------

  # A single element R vector is converted to a Python scalar. To overcome this
  # behavior R vectors are represented as lists explicitly.

  py_obj <- reticulate::py$convert_to_inchi(as.list(mols))
  inchi <- vapply(reticulate::py_to_r(py_obj), function(a1) {
    if(is.null(a1)) {
      return(NA_character_)
    } else {
      return(a1)
    }
  }, character(1L), USE.NAMES = FALSE)
  return(inchi)
}



#==============================================================================#
#' Convert molecular representations to InChIKey strings
#'
#' @description
#'   Convert molecular representations to InChIKey strings.
#'
#'   The input must be a character vector of either SMILES or InChI identifiers;
#'   mixing the two in the same vector is not allowed.
#'
#' @param mols
#'   A character vector. Molecular representation, this can be SMILES or InChI
#'   strings. All elements must be of the same type.
#'
#' @return
#'   A character vector. InChIKey strings. Elements that cannot be converted are
#'   returned as \code{NA}.
#'
#' @examples
#'   # Convert a vector of InChI to InChIKey identifiers
#'   inchi <- c("InChI=1S/C2H6/c1-2/h1-2H3",
#'              "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H")
#'   ConvertToInchikey(inchi)
#'   #> "OTMSDBZUPAUEDD-UHFFFAOYSA-N"
#'   #> "UHOVQNZJYSORNB-UHFFFAOYSA-N"
#'
#'   # Convert a vector of SMILES to InChIKey identifiers
#'   smiles <- c("CC", "c1ccccc1")
#'   ConvertToInchikey(smiles)
#'   #> "OTMSDBZUPAUEDD-UHFFFAOYSA-N"
#'   #> "UHOVQNZJYSORNB-UHFFFAOYSA-N"
#'
#'
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
ConvertToInchikey <- function(mols) {

  #--[ Pre-processing ]---------------------------------------------------------

  .EnsurePythonReady()
  mol_input_type <- .InferMolInputType(mols[1L])


  #--[ Check input arguments ]--------------------------------------------------

  if (is.na(mol_input_type)) {
    stop("'mols' must be a vector of SMILES or InChI identifiers.")
  }


  #--[ Convert InchI ]----------------------------------------------------------

  # A single element R vector is converted to a Python scalar. To overcome this
  # behavior R vectors are represented as lists explicitly.

  py_obj <- reticulate::py$convert_to_inchikey(as.list(mols))
  inchikey <- vapply(reticulate::py_to_r(py_obj), function(a1) {
    if(is.null(a1)) {
      return(NA_character_)
    } else {
      return(a1)
    }
  }, character(1L), USE.NAMES = FALSE)
  return(inchikey)
}


