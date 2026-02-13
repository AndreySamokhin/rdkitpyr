#==============================================================================#
#' @importFrom reticulate py_require
#==============================================================================#
.onLoad <- function(libname, pkgname) {

  the$py_requirements <- list(
    python = c(op = ">=", ver = "3.12"),
    packages = list(
      numpy = list(c(op = ">=",  ver = "2"),
                   c(op = "<",   ver = "3")),
      rdkit = list(c(op = ">=",  ver = "2025.9.1"))
    )
  )
  the$python_ready <- FALSE
  the$rdkit_available <- NA

  py_packages <- character(0L)
  for (py_pkg_name in names(the$py_requirements$packages)) {
    for (i in seq_along(the$py_requirements$packages[[py_pkg_name]])) {
      ver <- paste(the$py_requirements$packages[[py_pkg_name]][[i]],
                   collapse = "")
      py_packages[[length(py_packages) + 1L]] <- paste0(py_pkg_name, ver)
    }
  }
  reticulate::py_require(
    python_version = paste(the$py_requirements$python, collapse = ""),
    packages = py_packages
  )
  # cat(paste(py_packages, collapse = "\n"))
}


