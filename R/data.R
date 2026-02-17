#' Test compounds
#'
#' @description
#'   A small reference set of chemical compounds extracted from the PubChem
#'   database. Despite its limited size, the dataset covers a broad range of
#'   chemical properties, including:
#'   \itemize{
#'     \item{neutral, charged, and radical species;}
#'     \item{aromatic and aliphatic compounds;}
#'     \item{molecules containing heteroatoms;}
#'     \item{stereochemistry;}
#'     \item{species with disconnected fragments.}
#'   }
#'
#' @format
#'   A data frame.
#'   \describe{
#'     \item{\code{pubchem_cid}}{PubChem compound identifier (CID).}
#'     \item{\code{name}}{Compound name (as reported in the PubChem database).}
#'     \item{\code{smiles}}{SMILES identifier.}
#'     \item{\code{inchi}}{InChI identifier.}
#'     \item{\code{inchikey}}{InChIKey identifier.}
#'     \item{\code{n_labeled_atoms}}{Number of isotopically labeled atoms.}
#'     \item{\code{charge}}{The total charge of a molecule.}
#'     \item{\code{exact_mass}}{Exact monoisotopic mass.}
#'     \item{\code{molecular_weight}}{Average molecular weight.}
#'     \item{\code{formula}}{Molecular formula.}
#'     \item{\code{formula_isotopes}}{Molecular formula with explicit isotopes.}
#'   }
"test_compounds"


