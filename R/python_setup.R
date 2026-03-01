#==============================================================================#
#' Get information about the currently used Python environment
#'
#' @description
#'   Get details about the Python interpreter, \code{numpy} and \code{rdkit}
#'   packages.
#'
#' @param verbose
#'   A logical value. If \code{TRUE} (default), function prints information
#'   about currently used Python environemt.
#'
#' @return
#'   Invisibly return a named list with the following components:
#'   \describe{
#'     \item{python_path}{Full path to the Python executable.}
#'     \item{python_version}{Version of the Python interpreter.}
#'     \item{numpy_version}{Installed version of the \code{numpy} package.}
#'     \item{rdkit_version}{Installed version of the \code{rdkit} package.}
#'     \item{forced_by}{If applicable, indicates whether the Python
#'     interpreter was forced via \code{RETICULATE_PYTHON}, \code{use_*()}, or
#'     \code{py_require()}.}
#'   }
#'
#' @examplesIf rdkitpyr:::.IsRdkitAvailable(initialize = FALSE)
#'   # Print information about the Python environment
#'   GetPythonInfo()
#'
#'   # Access programmatically
#'   py_env <- GetPythonInfo(verbose = FALSE)
#'   py_env$python_version
#'   py_env$rdkit_version
#'
#' @export
#==============================================================================#
GetPythonInfo <- function(verbose = TRUE) {
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a logical value.")
  }
  .EnsurePythonReady()
  if (verbose) {
    forced_by <- ""
    if (!is.null(the$py_info$forced_by)) {
      forced_by <- paste0("\n\nPython version is forced by ",
                          the$py_info$forced_by)
    }
    message("Python path:    ", the$py_info$python_path, "\n",
            "Python version: ", the$py_info$python_version, "\n",
            "numpy version:  ", the$py_info$numpy_version, "\n",
            "rdkit version:  ", the$py_info$rdkit_version,
            forced_by)
  }
  return(invisible(the$py_info))
}



#==============================================================================#
#' Check whether RDKit is available
#'
#' @description
#'   Used in examples and tests to determine whether the Python module
#'   \code{rdkit} is available via \pkg{reticulate}.
#'
#'   This function is exported to support examples and tests. It is not part of
#'   the stable user-facing API and may change without notice.
#'
#'   The result is cached for the duration of the R session to avoid repeated
#'   calls to \code{reticulate::py_module_available()}.
#'
#' @param initialize
#'   A logical value. If \code{FALSE} and Python is not initialized, RDKit is
#'   assumed to be unavailable and the corresponding tests and examples are
#'   skipped.
#'
#' @return
#'   A logical value indicating whether the \code{rdkit} Python module is
#'   available.
#'
#' @importFrom reticulate py_available
#' @importFrom reticulate py_module_available
#'
#' @export
#==============================================================================#
.IsRdkitAvailable <- function(initialize = TRUE) {
  if (!is.na(the$rdkit_available)) {
    return(the$rdkit_available)
  }
  if (!initialize && !reticulate::py_available(initialize = FALSE)) {
    return(FALSE)
  }
  the$rdkit_available <- reticulate::py_module_available("rdkit")
  return(the$rdkit_available)
}



#==============================================================================#
#' Compare a version string with a required version
#' @param cur_ver A string. Current version.
#' @param op A string. Operator (e.g., \code{"=="}, \code{">="}, etc.).
#' @param req_ver A string. Required version.
#' @importFrom utils compareVersion
#' @noRd
#==============================================================================#
.VerifyVersion <- function(cur_ver, op, req_ver) {
  res <- utils::compareVersion(cur_ver, req_ver)
  out <- switch(op,
                "==" = (res == 0),
                "!=" = (res != 0),
                ">"  = (res >  0),
                ">=" = (res >= 0),
                "<"  = (res <  0),
                "<=" = (res <= 0),
                stop("Incorrect operator is used to compare two versions."))
  return(out)
}



#==============================================================================#
#' Check that Python meets the required version
#' @importFrom reticulate py_available
#' @importFrom reticulate py_config
#' @noRd
#==============================================================================#
.VerifyPythonVersion <- function() {
  if (!reticulate::py_available()) {
    stop("If this function is called, Python must already be initialized.")
  }
  py_conf <- reticulate::py_config()
  the$py_info$python_version <- as.character(py_conf$version)
  op <- the$py_requirements$python[["op"]]
  req_ver <- the$py_requirements$python[["ver"]]
  is_version_valid <- .VerifyVersion(the$py_info$python_version, op, req_ver)
  if (!is_version_valid) {
    warning("The following requirement is not fulfilled: ", "Python",
            paste(the$py_requirements$python, collapse = ""), call. = FALSE)
  }
  return(invisible(NULL))
}



#==============================================================================#
#' Check that Python packages are installed and meet version constraints
#' @importFrom reticulate py_available
#' @importFrom reticulate import
#' @noRd
#==============================================================================#
.VerifyPythonPackageVersions <- function() {
  if (!reticulate::py_available()) {
    stop("If this function is called, Python must already be initialized.")
  }
  importlib <- reticulate::import("importlib.metadata")
  for (py_pkg_name in names(the$py_requirements$packages)) {
    # If package is not installed, 'importlib$version(py_pkg_name)' returns the
    # following error: "No package metadata was found for <PKG_NAME>"
    py_pkg_ver <- importlib$version(py_pkg_name)
    if (!is.character(py_pkg_ver) || length(py_pkg_ver) != 1L) {
      the$py_info[[paste0(py_pkg_name, "_version")]] <- NA_character_
      warning("Unable to determine version of Python package '",
              py_pkg_name, "'.", call. = FALSE)
    } else {
      the$py_info[[paste0(py_pkg_name, "_version")]] <- py_pkg_ver
      for (i in seq_along(the$py_requirements$packages[[py_pkg_name]])) {
        req <- the$py_requirements$packages[[py_pkg_name]][[i]]
        if (!.VerifyVersion(py_pkg_ver, req[["op"]], req[["ver"]])) {
          warning("The following requirement is not fulfilled: ",
                  paste0(py_pkg_name, req[["op"]], req[["ver"]]), call. = FALSE)
        }
      }
    }
  }
  return(invisible(NULL))
}



#==============================================================================#
#' Load Python module
#' @importFrom reticulate py_run_string
#' @importFrom reticulate import
#' @noRd
#==============================================================================#
.LoadPythonModule <- function() {
  # The 'insert(0, path)' has higher prioriy over 'append(path)', so it prevents
  # conflicts and ensures that the module will be loaded.
  py_module_path <- system.file("python", package = "rdkitpyr")
  py_script <- sprintf("import sys; sys.path.insert(0, '%s')", py_module_path)
  reticulate::py_run_string(py_script)
  the$py_module <- reticulate::import("rdkitpyr", convert = FALSE)
  return(invisible(NULL))
}



#==============================================================================#
#' Initialize and validate Python environment for package use
#' @importFrom reticulate py_config
#' @noRd
#==============================================================================#
.EnsurePythonReady <- function() {

  # Requirements set with 'py_require()' are not taken into account, when Python
  # is forced manually, for example, by RETICULATE_PYTHON_ENV or 'use_*()'. In
  # such cases manual verification based on '.VerifyPythonVersion()' and
  # '.VerifyPythonPackageVersions()' is used.

  if (isTRUE(the$py_ready)) {
    return(invisible(NULL))
  }
  config <- reticulate::py_config() # explicit Python initialization
  if (!reticulate::py_available()) {
    stop("Python initialization failed.", call. = FALSE)
  }
  the$py_info$python_path <- config$python
  if (!is.null(config$forced)) {
    the$py_info$forced_by <- config$forced
  }
  .VerifyPythonVersion()
  .VerifyPythonPackageVersions()
  .LoadPythonModule()
  the$py_ready <- TRUE
  return(invisible(NULL))
}


