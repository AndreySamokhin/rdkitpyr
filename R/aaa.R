the <- new.env(parent = emptyenv())
the$py_requirements <- list(python = NULL,
                            packages = NULL)
the$py_info <- list(python_path    = NULL,
                    python_version = NULL,
                    numpy_version  = NULL,
                    rdkit_version  = NULL,
                    forced_by      = NULL)
the$py_ready <- FALSE

