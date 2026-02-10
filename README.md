# rdkitpyr

The `rdkitpyr` package provides a high-level interface for chemoinformatics
tasks in R. The RDKit library serves as a back end (accessed via its Python API
through the `reticulate` package). Rather than acting as a thin wrapper around
RDKit, `rdkitpyr` exposes R-friendly, task-oriented functions without mirroring
RDKit's native API structure. The package motto is *"Do chemistry tasks from R,
not use RDKit in R"*.

The name `rdkitpyr` (pronounced `/ar-dee-kit-pier/`) reflects both its nested
structure (from R to Python to the RDKit C++ library) and its role as an
R-side peer to RDKit, not a direct wrapper.

A key design concept of the package is one Python call per task. Custom Python
helper functions implement this approach. They reduce overhead, and handle
exceptions consistently.


## Installation

``` r
# Install 'rdkitpyr' from GitHub:
library(devtools)
install_github("https://github.com/AndreySamokhin/rdkitpyr")
```


## Quick Start

The `rdkitpyr` package supports both novice and advanced users.

* Advanced users can manually set up a Python environment (`Python>=3.12`,
`numpy>= 2,<3`, and `rdkit>=2025.9.1`) and select it in R using one of the
methods recommended by `reticulate` (e.g., setting the `RETICULATE_PYTHON_ENV`
environment variable or using `reticulate::use_*()` functions). This provides
full control over package versions and the Python environment.
* Novice users do not need to install or configure Python manually. This is
handled automatically by `reticulate` using the manifest created with
`reticulate::py_require()`. Installation starts automatically when a function is
first called. It may take a few minutes, and requires about 200 MB of disk
space.

Here is a simple example:

```r
library(rdkitpyr)

rdkitpyr::ConvertToInchi(c("CC", "CCC"))
#> "InChI=1S/C2H6/c1-2/h1-2H3"
#> "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3"
```


