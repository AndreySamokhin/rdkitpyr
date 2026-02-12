#==============================================================================#
#' Convert molecular representations to SMILES strings
#'
#' @description
#'   Convert molecular representations to SMILES strings.
#'
#'   The input must be a character vector of SMILES or InChI strings,
#'   or a list of RDKit Mol objects.
#'
#'   When SMILES strings are provided as input, their representation can
#'   be transformed, for example, to canonical SMILES or with stereochemistry
#'   omitted.
#'
#' @inheritParams ParseMolecules
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
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
                            canonical = TRUE,
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

  # 'verbose'
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a logical value.")
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
#' @inheritParams ParseMolecules
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
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
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
ConvertToInchi <- function(mols,
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


  #--[ Convert Smiles ]---------------------------------------------------------

  # A single element R vector is converted to a Python scalar. To overcome this
  # behavior R vectors are represented as lists explicitly.

  py_obj <- reticulate::py$convert_to_inchi(as.list(mols), verbose = verbose)
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
#'   The input must be a character vector of SMILES or InChI strings,
#'   or a list of RDKit Mol objects.
#'
#' @inheritParams ParseMolecules
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
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
#' @importFrom reticulate py
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
ConvertToInchikey <- function(mols,
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



#==============================================================================#
#' Parse SMILES and InChI strings into RDKit Mol objects
#'
#' @description
#'   Parse SMILES and InChI strings into RDKit Mol objects
#'
#'   This function converts a character vector of molecular representations
#'   (SMILES or InChI) into a list of RDKit Mol objects (Python-backed pointers
#'   via \code{reticulate}). The resulting objects can be reused in subsequent
#'   operations without repeated conversion from SMILES or InChI. This is
#'   particularly useful when multiple chemoinformatics tasks are performed on
#'   the same set of molecules, improving efficiency by avoiding repeated
#'   parsing steps.
#'
#' @param mols
#'   A character vector of SMILES or InChI strings.
#' @param verbose
#'   A logical value. If \code{TRUE}, display warning and error messages emitted
#'   by the RDKit backend. Default is \code{FALSE}.
#'
#' @return
#'   A list of RDKit Mol objects.
#'
#' @examples
#'   # Convert a vector of SMILES to RDKit Mol objects
#'   mols <- ParseMolecules(c("CC", "CCC"))
#'   print(mols[[1L]])
#'   #> <rdkit.Chem.rdchem.Mol object at 0x000001CC4D60F4C0>
#'
#'   # Convert a list of RDKit Mol objects to InChI identifiers
#'   ConvertToInchi(mols)
#'   #> "InChI=1S/C2H6/c1-2/h1-2H3"
#'   #> "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3"
#'
#' @importFrom reticulate py
#'
#' @export
#==============================================================================#
ParseMolecules <- function(mols,
                           verbose = FALSE) {
  .EnsurePythonReady()
  if (!is.character(mols) || length(mols) == 0L) {
    stop("'mols' must be a character vector of SMILES or InChI strings.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a logical value.")
  }
  return(reticulate::py$convert_to_molecules(as.list(mols)))
}


