# rdkitpyr

The `rdkitpyr` package provides a high-level interface for cheminformatics
tasks in R. It leverages the RDKit library via its Python API, accessed through
the `reticulate` package. Unlike a thin wrapper around RDKit, `rdkitpyr` offers
R-friendly, task-oriented functions without mirroring RDKit's native API. The
package motto is *"Focus on chemistry tasks in R, without diving into RDKit
internals"*.

The name `rdkitpyr` (pronounced `/ar-dee-kit-peer/`) highlights its nested
structure (from R to Python to the RDKit C++ library) and its role as an R-side
peer to RDKit rather than a direct wrapper.

A key design principle is one Python call per task. Custom Python helper
functions implement this approach, reducing overhead and ensuring consistent
exception handling.


## Installation

``` r
# Install 'rdkitpyr' from CRAN:
install.packages("rdkitpyr")
```

``` r
# Install 'rdkitpyr' from GitHub:
library(devtools)
install_github("andreysamokhin/rdkitpyr")
```


## Quick Start

The `rdkitpyr` package supports both novice and advanced users.

* Advanced users can manually set up a Python environment (`Python>=3.12`,
`numpy>=2,<3`, and `rdkit>=2025.9.1`) and select it in R using one of the
methods recommended by `reticulate` (e.g., setting the `RETICULATE_PYTHON_ENV`
environment variable or using `reticulate::use_*()` functions). This gives
full control over package versions and the Python environment.
* Novice users do not need to install or configure Python manually. This is
handled automatically by `reticulate` using a manifest created with
`reticulate::py_require()`. Installation starts only when a function that
requires Python is called. It may take a few minutes, and requires about 200 MB
of disk space.

Here is a simple example:

```r
library(rdkitpyr)

rdkitpyr::ConvertToInchi(c("CC", "CCC"))
#> "InChI=1S/C2H6/c1-2/h1-2H3"
#> "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3"

rdkitpyr::CalculateExactMass(c("c1ccccc1", "CCCC"))
#> 78.04695 58.07825
```


