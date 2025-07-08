# Load required packages
library(tidyr)
library(dplyr)

# Function to convert CSV to array (keep your existing function)
csv_to_array <- function(file_path, has_bootstrap = FALSE) {
  # Read the CSV
  data <- read.csv(file_path)

  if (has_bootstrap) {
    # 4D array with bootstrap
    dims <- c(
      max(data$l),  # Bootstrap samples
      max(data$i),  # Outcomes
      max(data$j),  # Subgroups
      max(data$k)   # Comparisons
    )

    # Create empty array
    result <- array(NA, dim = dims)

    # Fill the array
    for (row in 1:nrow(data)) {
      l <- data$l[row]
      i <- data$i[row]
      j <- data$j[row]
      k <- data$k[row]
      result[l, i, j, k] <- data$value[row]
    }
  } else {
    # 3D array without bootstrap
    dims <- c(
      max(data$i),  # Outcomes
      max(data$j),  # Subgroups
      max(data$k)   # Comparisons
    )

    # Create empty array
    result <- array(NA, dim = dims)

    # Fill the array
    for (row in 1:nrow(data)) {
      i <- data$i[row]
      j <- data$j[row]
      k <- data$k[row]
      result[i, j, k] <- data$value[row]
    }
  }

  return(result)
}

# Load the Stata bootstrap data
abregact <- csv_to_array("/Users/elliotpaschal/Documents/mhtexp2r/Tests/abregact.csv")
abregboot <- csv_to_array("/Users/elliotpaschal/Documents/mhtexp2r/Tests/abregboot.csv", has_bootstrap = TRUE)
pact <- csv_to_array("/Users/elliotpaschal/Documents/mhtexp2r/Tests/pact.csv")
pboot <- csv_to_array("/Users/elliotpaschal/Documents/mhtexp2r/Tests/pboot.csv", has_bootstrap = TRUE)

# Create combo matrix
combo_mat <- matrix(c(0, 1, 0, 2, 1, 2), nrow = 3, byrow = TRUE)

# Create select array (assuming all hypotheses are selected)
num_outcomes <- dim(abregact)[1]
num_subgroups <- dim(abregact)[2]
num_comparisons <- dim(abregact)[3]
select <- array(1, dim = c(num_outcomes, num_subgroups, num_comparisons))

# Use the Stata p-values directly (no need to recalculate)
pvals <- pact  # Use Stata's calculated p-values directly

# Calculate all the threshold corrections using your updated functions
# Use Stata's pboot directly for alphasin calculation
alpha1 <- calculate_alphasin(pvals, pboot)

# Add this debug right before calling calculate_alpha_unified:
cat("Array dimensions check:\n")
cat("pvals dimensions:", dim(pvals), "\n")
cat("pboot dimensions:", dim(pboot), "\n")
# cat("results$boot dimensions:", dim(results$boot), "\n")
# cat("coefficients dimensions:", dim(results$coef), "\n")

# Use your unified threshold function
thresholds <- calculate_alpha_unified(
  pvals = pvals,
  coefficients = abregact,  # Use the loaded coefficients
  pboot = pboot,
  combo = combo_mat,
  alphasin = alpha1,
  select = select,
  transitivitycheck = TRUE
)

alpha2 <- thresholds$alphamul   # Standard thresholds
alpha3 <- thresholds$alphamulm  # Transitivity-aware thresholds

# Calculate Bonferroni and Holm corrections
num_hypotheses <- num_outcomes * num_subgroups * num_comparisons
alphasin_flat <- as.vector(alpha1)

# Bonferroni correction
bonf_flat <- pmin(1, alphasin_flat * num_hypotheses)

# Holm correction
o <- order(alphasin_flat)
alphasin_sorted <- alphasin_flat[o]
holm_adjust <- alphasin_sorted * rev(seq_len(num_hypotheses))
holm_adjust <- pmin(cummax(holm_adjust), 1)
holm_flat <- numeric(num_hypotheses)
holm_flat[o] <- holm_adjust

# Create results manually to match Stata ordering
# Create indices in the same order as Stata (outcome varies slowest)
indices <- expand.grid(
  comparison = 1:num_comparisons,
  subgroup = 1:num_subgroups,
  outcome = 1:num_outcomes
)

# Reorder to match your R convention (outcome varies slowest)
indices <- indices[order(indices$outcome, indices$subgroup, indices$comparison), ]

# Add treatment pairs
treatment_pairs <- t(sapply(indices$comparison, function(comp_idx) {
  c(combo_mat[comp_idx, 1], combo_mat[comp_idx, 2])
}))
colnames(treatment_pairs) <- c("t1", "t2")

# Create results data frame manually
results_output <- data.frame(
  outcome = indices$outcome,
  subgroup = indices$subgroup,
  comparison = indices$comparison,
  treatment_pairs,
  coefficient = as.vector(abregact),  # Using test stats as proxy
  test_stat = as.vector(abregact),
  Remark3_2 = as.vector(alpha1),
  Thm3_1 = as.vector(alpha2),
  Remark3_8 = as.vector(alpha3),
  Bonf = bonf_flat,
  Holm = holm_flat
)

cat("\nResults using Stata bootstrap data with corrected ordering:\n")
print(results_output, digits = 4)

# Compare with expected Stata output
cat("\nKey comparisons:\n")
cat("Remark3_2 (first few):", head(as.vector(alpha1)), "\n")
cat("Thm3_1 (first few):", head(as.vector(alpha2)), "\n")
cat("Remark3_8 (first few):", head(as.vector(alpha3)), "\n")
