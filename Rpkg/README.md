# ORION R Package

This R package provides data objects and minimal R infrastructure associated with the ORION project. 
It serves as a lightweight container for experiment metadata, processed datasets, and documentation 
used in downstream analyses. No modeling or computational routines are included.

## Installation

You can install the development version directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("eriklarsen4/ORION", subdir = "Rpkg")
```

## Contents

- `data/` — packaged datasets
- `vignettes/` — documentation and examples  

## Usage

Load the package:

```r
library(orion)
```

Then access packaged data:

```r
data("example_dataset")
```

## Project Structure

This R package is part of the larger ORION repository, which also contains:

- Python modeling code (`python/orion`)  
- Analysis scripts (`analysis/`)  
- Docker configurations (`docker/`)  
