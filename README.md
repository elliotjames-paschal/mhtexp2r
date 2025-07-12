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
devtools::install_github("elliotjames-paschal/mhtexp2r")
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

## Validation

The package includes a validation function to compare the R implementation against the original Stata `mhtexp2` command:

```r
# Basic validation with default settings
result <- validate_mhtexp2()

# The function will:
# - Auto-detect your Stata installation
# - Generate 1000 test observations with 2 outcomes, 2 subgroups, 3 treatments
# - Run 5000 bootstrap replications for precise p-value estimation
# - Compare R vs Stata results across all correction methods

# View the comparison table
View(result$comparison)

# Customized validation
result <- validate_mhtexp2(
  n_obs = 2000,                    # More observations
  bootstrap = 10000,               # More bootstrap replications  
  combo = "treatmentcontrol",      # Different comparison type
  seed = 42,                       # Reproducible results
  stata_path = "/custom/path/to/stata",  # Manual Stata path if auto-detection fails
  verbose = TRUE                   # Show detailed output
)

# Use your own data
result <- validate_mhtexp2(
  data_source = "path/to/your/data.csv",  # Must have columns: y1, y2, treat, x1, x2, subgroup
  bootstrap = 5000,
  stata_path = "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"  # Specify Stata location
)

# Check what Stata installation was found
find_stata()  # Returns path to detected Stata executable
```

**Requirements:**
- Stata installation (any recent version)
- Data must include columns: `y1`, `y2`, `treat`, `x1`, `x2`, `subgroup`

**Note:** If Stata auto-detection fails, manually specify the path with `stata_path = "/path/to/your/stata"`.

The validation compares coefficients and all p-value correction methods (Remark 3.2, Theorem 3.1, Remark 3.8) between R and Stata implementations.

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

