# mhtexp2r

**Multiple Hypothesis Testing for Experimental Economics**

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`mhtexp2r` implements multiple hypothesis testing procedures for experimental economics with exact replication of `mhtexp2` functionality from List, Shaikh, and Vayalinkal (2019). The package provides family-wise error rate control for treatment effect testing with bootstrap inference and transitivity corrections.

## Installation

You can install the package directly from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install mhtexp2r from GitHub
devtools::install_github("yourusername/mhtexp2r")
```

## Usage

### Basic Example

```r
library(mhtexp2r)

# Load your experimental data
# Assume: df with columns y1, y2 (outcomes), treat (treatment), x1, x2 (controls), subgroup

# Run multiple hypothesis testing
results <- mhtexp2_r(
  Y = df[, c("y1", "y2")],           # Matrix of outcomes
  treatment = df$treat,               # Treatment assignment vector
  controls = df[, c("x1", "x2")],    # Control variables (optional)
  subgroup = df$subgroup,             # Subgroup identifiers (optional)
  combo = "pairwise",                 # Comparison type: "pairwise" or "treatmentcontrol"
  bootstrap = 3000,                   # Number of bootstrap replications
  studentized = TRUE,                 # Whether to studentize test statistics
  transitivity_check = TRUE           # Apply transitivity corrections
)

# View results
head(results$output)
```

### Output Description

The main output is a data frame with the following columns:

- **`outcome`**: Outcome variable index
- **`subgroup`**: Subgroup index  
- **`comparison`**: Treatment comparison index
- **`t1`, `t2`**: Treatment pair being compared
- **`coefficient`**: Treatment effect estimate
- **`test_stat`**: Test statistic
- **`Remark3_2`**: Single hypothesis adjusted p-values
- **`Thm3_1`**: Multiple hypothesis adjusted p-values  
- **`Remark3_8`**: Transitivity-corrected p-values
- **`Bonf`**: Bonferroni correction
- **`Holm`**: Holm correction

## Testing

The package includes validation tests that compare R and Stata implementations:

### Validation Script

Use `validate_package.R` to run both implementations side-by-side:

```r
# Load and run validation
result <- validate_mhtexp2()

# View detailed comparison
View(result$comparison)
```

### Additional Test Files

- **`export_bootstrap.do`** - Modified Stata script that exports intermediate bootstrap values
- **`statatest.R`** - Functions to import and compare detailed Stata exports with R calculations
- **Standard tests** - Located in `tests/testthat/`, run with `devtools::test()`

## Requirements

- R >= 4.0.0
- No additional dependencies for core functionality
- Optional: `dplyr`, `tidyr`, `testthat` for testing

## References

List, J. A., Shaikh, A. M., & Vayalinkal, J. P. (2019). Multiple Testing with Covariate Adjustment in Experimental Economics. *Journal of Econometrics*, forthcoming.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

elliotjpaschal@gmail.com

