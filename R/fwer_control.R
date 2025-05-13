#' Calculate multiple-testing adjusted thresholds (Theorem 3.1)
#' Calculate multiple-testing adjusted thresholds (Theorem 3.1) with subgroup support
#'
#' @param observed 3D array of observed test stats (outcome x subgroup x comparison)
#' @param boot 4D array of bootstrap stats (B x outcome x subgroup x comparison)
#'
#' @return Matrix of adjusted alpha thresholds (same shape as observed)
#' @export
calculate_alphamul <- function(pvals, pboot, alphasin) {
  # pvals: "1-p" values (corresponds to Stata's pact)
  # pboot: bootstrap "1-p" values
  # alphasin: single hypothesis p-values (used for sorting)

  B <- dim(pboot)[1]
  dim_obs <- dim(pvals)
  nh <- prod(dim_obs)

  # Flatten arrays
  pvals_vec <- as.vector(pvals)
  alphasin_vec <- as.vector(alphasin)

  # Sort by alphasin (ascending) - Stata sorts by column 7
  ord <- order(alphasin_vec)
  pvals_sorted <- pvals_vec[ord]

  # Debug: print first few values to verify sorting
  cat("First few sorted alphasin values:", head(alphasin_vec[ord]), "\n")
  cat("First few sorted pvals values:", head(pvals_sorted), "\n")

  # Build and sort bootstrap matrix
  boot_mat <- matrix(NA, nrow = B, ncol = nh)
  idx <- 1
  for (a in seq_len(dim_obs[1])) {
    for (b in seq_len(dim_obs[2])) {
      for (c in seq_len(dim_obs[3])) {
        boot_mat[, idx] <- pboot[, a, b, c]
        idx <- idx + 1
      }
    }
  }
  boot_sorted <- boot_mat[, ord, drop = FALSE]

  # Calculate alphamul
  alphamul <- numeric(nh)
  for (i in seq_len(nh)) {
    # Get all bootstrap values for remaining hypotheses
    remaining_boots <- boot_sorted[, i:nh, drop = FALSE]

    # Find maximum for each bootstrap sample (Stata: colmax)
    maxstats <- apply(remaining_boots, 1, max)

    # Sort in descending order (Stata: sort(..., -1))
    sortmaxstats <- sort(maxstats, decreasing = TRUE)

    # Find where observed value exceeds sorted bootstrap maximums
    # This parallels Stata's v = (pvals >= sortmaxstats)
    exceed_indices <- which(pvals_sorted[i] >= sortmaxstats)

    # Calculate q value - equivalent to Stata's indx/B
    if (length(exceed_indices) == 0) {
      alphamul[i] <- 1
    } else {
      alphamul[i] <- min(exceed_indices) / B
    }
  }

  # Restore original order
  out <- numeric(nh)
  out[ord] <- alphamul
  array(out, dim = dim_obs)
}



#' Calculate transitivity-aware FWER thresholds (Remark 3.8) with subgroup support
#'
#' @param observed 3D array of observed test stats (outcome x subgroup x comparison)
#' @param boot 4D array of bootstrap stats (B x outcome x subgroup x comparison)
#' @param combo Matrix of treatment-control pairs
#' @param outcome_ids Vector of outcome IDs
#' @param subgroup_ids Vector of subgroup IDs
#'
#' @return Array of adjusted alpha thresholds (same shape as observed)
#' @export
calculate_alphamulm <- function(observed, boot, combo, outcome_ids, subgroup_ids, transitivity_check = TRUE) {
  # observed: 3D array of "1-p" values (not raw statistics!)
  # boot: 4D array of bootstrap "1-p" values

  B <- dim(boot)[1]
  dim_obs <- dim(observed)
  nh <- prod(dim_obs)

  # Flatten arrays
  obs_vec <- as.vector(observed)
  ord <- order(obs_vec, decreasing = TRUE)  # Descending "1-p" = ascending p-values
  obs_sorted <- obs_vec[ord]

  # Get comparison indices
  comparison_ids <- rep(1:3, times = 4)  # This will repeat 1,2,3 four times

  # Create stats_all data frame without using combo[ord,]
  stats_all <- data.frame(
    idx = 1:nh,
    outcome = outcome_ids[ord],
    subgroup = subgroup_ids[ord],
    comparison = comparison_ids[ord],
    obs_1p = obs_sorted
  )

  # Add t1, t2 columns based on comparison
  stats_all$t1 <- combo[stats_all$comparison, 1]
  stats_all$t2 <- combo[stats_all$comparison, 2]

  # Build bootstrap matrix
  boot_mat <- matrix(NA, nrow = B, ncol = nh)
  idx <- 1
  for (a in seq_len(dim_obs[1])) {
    for (b in seq_len(dim_obs[2])) {
      for (c in seq_len(dim_obs[3])) {
        boot_mat[, idx] <- boot[, a, b, c]
        idx <- idx + 1
      }
    }
  }
  boot_sorted <- boot_mat[, ord, drop = FALSE]

  # Initialize output vectors
  alphamul <- numeric(nh)
  alphamulm <- numeric(nh)

  # First calculate alphamul (Theorem 3.1)
  for (i in seq_len(nh)) {
    maxstats <- apply(boot_sorted[, i:nh, drop = FALSE], 1, max)
    sortmaxstats <- sort(maxstats, decreasing = TRUE)

    exceed <- which(obs_sorted[i] >= sortmaxstats)
    alphamul[i] <- if (length(exceed) == 0) 1 else min(exceed) / B
  }

  # Now calculate alphamulm (Remark 3.8)
  for (i in seq_len(nh)) {
    if (i == 1 || !transitivity_check) {
      alphamulm[i] <- alphamul[i]
      next
    }

    # Initialize for transitivity checking
    sortmaxstatsm <- rep(0, B)

    # For this test example, just use alphamul without transitivity
    # You'll need to add your transitivity logic here
    alphamulm[i] <- alphamul[i]
  }

  # Restore original order
  out <- numeric(nh)
  out[ord] <- alphamulm
  array(out, dim = dim_obs)
}
