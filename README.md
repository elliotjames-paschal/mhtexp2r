
---

# mhtexp2r

An R implementation of the `mhtexp2` multiple hypothesis testing procedure for experimental economics, replicating the Stata version developed in John List’s lab.

## Installation

### Option 1: Install from GitHub (recommended for users)

```r
install.packages("devtools")
devtools::install_github("elliotjames-paschal/mhtexp2r")
library(mhtexp2r)
```

### Option 2: Load locally (for developers)

If you cloned or downloaded the project folder directly:

```r
install.packages("devtools")
devtools::load_all(".")
```

## Usage Example

```r
df <- read.csv("data/testdata.csv")

results <- mhtexp2_r(
  Y = df[, c("y1", "y2")],
  treatment = df$treat,
  controls = df[, c("x1", "x2")],
  subgroup = df$subgroup,
  combo = "pairwise",
  bootstrap = 20,
  studentized = TRUE,
  transitivity_check = TRUE
)

head(results)
```

##  Project Structure

```
R/
├── mhtexp2_r.R           # Main interface
├── bootstrap_engine.R    # Bootstrap resampling and subgrouping
├── runreg.R              # Regression logic
├── input_builders.R      # Combo, select, subgroup parsers
├── pval_helpers.R        # P-value and alpha calculations
├── transitivity_check.R  # Optional Remark 3.8 logic
├── fwer_control.R        # Bonferroni & Holm methods
├── build_output.R        # Output shaping
```
