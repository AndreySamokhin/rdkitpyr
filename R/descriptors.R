#==============================================================================#
#' Calculate all RDKit molecular descriptors
#'
#' @description
#'   Calculate all molecular descriptors available in RDKit for a set of
#'   molecules.
#'
#'   The descriptors are calculated using the \code{CalcMolDescriptors()}
#'   function from the \code{Descriptors} module in RDKit.
#'
#'   The set of returned descriptors may depend on the installed RDKit version.
#'
#'   Each molecule is represented by a full set of descriptor values returned
#'   as a data frame. Invalid molecules are represented by rows containing
#'   \code{NA} values. Row order is preserved so that the output aligns with
#'   the input.
#'
#' @param mols
#'   A character vector of SMILES or InChI strings, or a list of RDKit
#'   Mol objects.
#' @inheritParams ParseMolecules
#'
#' @return
#'   A data frame with one row per molecule and one column per descriptor.
#'   Elements corresponding to invalid molecules are returned as \code{NA}.
#'   Original RDKit descriptor names are used.
#'
#'   Additionally, the \code{"valid"} attribute is attached to indicate which
#'   molecules were successfully processed.
#'
#' @examplesIf rdkitpyr:::.IsRdkitAvailable(initialize = FALSE)
#'   # Calculate all RDKit descriptors
#'   smiles <- c("CCO", "c1ccccc1", "invalid_molecule")
#'   desc <- CalculateAllDescriptors(smiles)
#'
#'   # Inspect names of first three descriptors
#'   names(desc)[1:3]
#'   #> "MaxAbsEStateIndex" "MaxEStateIndex" "MinAbsEStateIndex"
#'
#'   # Display Molecular weight, LogP, TPSA, and molar refractivity
#'   desc[c("MolWt", "MolLogP", "TPSA", "MolMR")]
#'   #>    MolWt MolLogP  TPSA   MolMR
#'   #> 1 46.069 -0.0014 20.23 12.7598
#'   #> 2 78.114  1.6866  0.00 26.4420
#'   #> 3     NA      NA    NA      NA
#'
#'   # Check which molecules were successfully processed
#'   attr(desc, "valid")
#'   #> TRUE TRUE FALSE
#'
#' @importFrom reticulate py_to_r
#'
#' @export
#==============================================================================#
CalculateAllDescriptors <- function(mols,
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

  py_obj <- the$py_module$calculate_all_descriptors(
    as.list(mols),
    verbose = verbose
  )
  temp <- reticulate::py_to_r(py_obj)
  n_descriptors <- length(temp[[1L]])
  descriptors <- as.data.frame(
    t(vapply(temp, function(a1) {
      num_vec <- unlist(a1)
      num_vec[is.nan(num_vec)] <- NA_real_
      return(num_vec)
    }, numeric(n_descriptors)))
  )
  attr(descriptors, "valid") <- !is.na(descriptors[, 1L])
  return(descriptors)
}



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
#' @examplesIf rdkitpyr:::.IsRdkitAvailable(initialize = FALSE)
#'   # Calculate exact mass for a set of molecules
#'   smiles <- c("CCO", "c1ccccc1", "invalid_molecule")
#'   CalculateExactMass(smiles)
#'   #> 46.04186 78.04695 NA
#'
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

  py_obj <- the$py_module$calculate_exact_mass(
    as.list(mols),
    verbose = verbose
  )
  out <- reticulate::py_to_r(py_obj)
  out[is.nan(out)] <- NA_real_
  return(out)
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
#' @examplesIf rdkitpyr:::.IsRdkitAvailable(initialize = FALSE)
#'   # Calculate average molecular weight for a set of molecules
#'   smiles <- c("CCO", "c1ccccc1", "invalid_molecule")
#'   CalculateMolecularWeight(smiles)
#'   #> 46.069 78.114 NA
#'
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

  py_obj <- the$py_module$calculate_molecular_weight(
    as.list(mols),
    verbose = verbose
  )
  out <- reticulate::py_to_r(py_obj)
  out[is.nan(out)] <- NA_real_
  return(out)
}


