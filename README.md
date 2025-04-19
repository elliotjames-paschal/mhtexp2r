
---

# mhtexp2r

An R implementation of the `mhtexp2` multiple hypothesis testing procedure for experimental economics, replicating the Stata version developed in John List’s lab.
## 🧠 Package Pipeline (Execution Flow)

The `mhtexp2r()` function performs subgroup-aware multiple hypothesis testing with bootstrap-based inference. It reproduces the behavior of the original Stata program line-for-line. Below is a detailed walkthrough of how the pipeline works, with each function labeled by the file it comes from.

---

### Step 1: **Input validation and setup**
**File:** `mhtexp2_r.R`

- Accepts user inputs: outcome matrix `Y`, treatment vector, optional controls, subgroup IDs, and testing parameters.
- Converts all inputs to matrices.
- Builds the treatment comparison matrix using `build_combo()`.
- Creates the inclusion `select` array for optional filtering of hypotheses (`only`, `exclude`).
- Generates or accepts a bootstrap ID matrix.

---

### Step 2: **Construct treatment comparisons**
**File:** `input_builders.R`

- `build_combo(groups, method = "pairwise" | "treatmentcontrol")`
  - If `"pairwise"`, builds all 2-way comparisons across groups.
  - If `"treatmentcontrol"`, compares each treatment to control (group 0).
  - Returns an n-row matrix of pairs.

---

### Step 3: **Bootstrap estimation loop**
**File:** `bootstrap_engine.R`

- `bootstrap_runreg()`:
  - Loops over outcomes, subgroups, and treatment pairs.
  - For each, estimates treatment effects using the `runreg()` function.
  - Repeats B times using bootstrap resampling (`idbootmat`).
  - Stores observed and simulated test statistics in a 4D array.
  - Applies `select[i, j, k]` logic to filter comparisons (optional).

**File:** `runreg.R`

- `runreg()`:
  - Computes adjusted ATE:  
    $begin:math:display$
    \\hat{ATE} = \\hat{\\beta}_1^{\\text{treated}} - \\hat{\\beta}_1^{\\text{control}} + \\bar{X}(b_1 - b_0)
    $end:math:display$
  - Estimates standard error using both heteroskedasticity and covariate adjustments.
  - Returns ATE and SE.

---

### Step 4: **P-value and alpha threshold calculations**
**File:** `pval_helpers.R`

- `calculate_pvals()`:
  - Computes the empirical p-value for each observed test stat based on bootstrap distribution.

- `calculate_alphasin()`:
  - Implements Remark 3.2: single-hypothesis rejection thresholds.

- `calculate_alphamul()`:
  - Implements Theorem 3.1: step-down multiple hypothesis testing thresholds (no transitivity).

---

### Step 5: **Transitivity correction (optional)**
**File:** `transitivity_check.R`

- `calculate_alphamulm()`:
  - Implements Remark 3.8 from Seidel and Xu (2020), ensuring no logical contradictions under the null when transitivity of rejected hypotheses is violated.
  - Used only if `transitivity_check = TRUE`.

---

### Step 6: **Output construction**
**File:** `build_output.R`

- `build_output()`:
  - Flattens and combines all arrays.
  - Outputs a `data.frame` with one row per hypothesis, including:
    - outcome index
    - subgroup
    - treatment and control group
    - test statistic
    - p-value
    - rejection thresholds from:
      - single-hypothesis testing
      - multiple-hypothesis testing (step-down)
      - transitivity-corrected step-down
    - Bonferroni and Holm p-values (from `fwer_control.R`)

---

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
