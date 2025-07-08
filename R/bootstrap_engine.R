#' Bootstrap adjusted regression with subgroup support
#' @param Y Matrix of outcomes (n x outcomes)
#' @param D Matrix (n x 1) of treatment assignments
#' @param DX Matrix (n x k) of treatment + controls
#' @param combo Matrix of treatment-control pairs
#' @param B Integer, number of bootstrap samples
#' @param subgroup Vector of subgroup IDs
#' @param studentized Whether to divide by SE
#' @param idbootmat Bootstrap resample matrix (n x B)
#' @param select 3D array specifying which outcome-subgroup-combo tests to run
#'
#' @return List with `stat` (test statistics), `coef` (ATEs), and `boot` arrays
#' @export
bootstrap_runreg <- function(Y, D, DX, combo, B,
                             subgroup = NULL,
                             studentized = TRUE,
                             idbootmat = NULL,
                             select = NULL) {

  # Quadvariance helper function
  quadvariance <- function(x) {
    if (is.null(x) || length(x) <= 1) return(0)

    if (is.matrix(x)) {
      # For matrices, apply quadvariance to each column
      n <- nrow(x)
      result <- matrix(0, ncol(x), ncol(x))
      means <- colMeans(x)

      # Calculate covariance matrix with n divisor
      for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
          result[i, j] <- sum((x[,i] - means[i]) * (x[,j] - means[j])) / n
        }
      }
      return(result)
    } else {
      # For vectors
      n <- length(x)
      mean_x <- mean(x)
      return(sum((x - mean_x)^2) / n)
    }
  }


  n <- nrow(Y)
  num_outcomes <- ncol(Y)
  num_combos <- nrow(combo)

  if (is.null(subgroup)) subgroup <- rep(1, n)
  subgroups <- sort(unique(subgroup))
  num_subgroups <- length(subgroups)

  observed_stat <- array(NA, dim = c(num_outcomes, num_subgroups, num_combos))
  observed_coef <- array(NA, dim = c(num_outcomes, num_subgroups, num_combos))
  boot_stats <- array(NA, dim = c(B, num_outcomes, num_subgroups, num_combos))

  # First calculate observed statistics
  for (i in seq_len(num_outcomes)) {
    for (s in seq_along(subgroups)) {
      sg <- subgroups[s]
      idx_sg <- which(subgroup == sg)

      Y_sg <- Y[idx_sg, , drop = FALSE]
      D_sg <- D[idx_sg, , drop = FALSE]
      DX_sg <- DX[idx_sg, , drop = FALSE]
      yi <- Y_sg[, i]

      for (j in seq_len(num_combos)) {
        if (!is.null(select) && select[i, s, j] != 1) next

        t1 <- combo[j, 1]
        t2 <- combo[j, 2]
        keep <- D_sg[,1] %in% c(t1, t2)


        if (sum(keep) < 2) next

        treat <- as.integer(D_sg[keep, 1] == t2)
        cur_X <- DX_sg[keep, , drop = FALSE]
        cur_X[,1] <- treat
        cur_y <- yi[keep]

        # Calculate controls for FULL subgroup (like Stata)
        if (ncol(DX_sg) > 1) {
          controls_full_sg <- DX[idx_sg, -1, drop = FALSE]
          Xbar_sg <- colMeans(controls_full_sg)
          n_sg <- nrow(controls_full_sg)
          Xvar_sg <- quadvariance(controls_full_sg)
        } else {
          Xbar_sg <- 0
          Xvar_sg <- 0
        }

        # Use original dataset for treatment probabilities (like Stata)
        pi_2 <- sum(subgroup == sg & D[,1] == t2) / n
        pi_1 <- sum(subgroup == sg & D[,1] == t1) / n
        pi_z <- sum(subgroup == sg) / n

        est <- runreg(cur_X, cur_y, Xbar_sg, Xvar_sg, pi_2, pi_1, pi_z)

        stat <- abs(est["ATE"]) / (if (studentized) est["SE"] else 1)

        observed_coef[i,s,j] <- est["ATE"]
        observed_stat[i,s,j] <- stat
      }
    }
  }

  # Now bootstrap
  for (b in seq_len(B)) {

    # First resample the dataset
    boot_idx <- idbootmat[,b]
    Y_boot <- Y[boot_idx, , drop = FALSE]
    D_boot <- D[boot_idx, , drop = FALSE]
    DX_boot <- DX[boot_idx, , drop = FALSE]
    sub_boot <- subgroup[boot_idx]

    # Process each outcome, subgroup, and comparison within this resampled dataset
    for (i in seq_len(num_outcomes)) {
      for (s in seq_along(subgroups)) {
        sg <- subgroups[s]
        idx_sg <- which(sub_boot == sg)

        if (length(idx_sg) < 2) next

        # Calculate controls for FULL bootstrap subgroup
        if (ncol(DX_boot) > 1) {
          controls_boot_sg <- DX_boot[idx_sg, -1, drop = FALSE]
          Xbar_b <- colMeans(controls_boot_sg)
          Xvar_b <- quadvariance(controls_boot_sg)

        } else {
          Xbar_b <- 0
          Xvar_b <- 0
        }

        # Subset the bootstrap sample for this outcome and subgroup
        Y_sg <- Y_boot[idx_sg, , drop = FALSE]
        D_sg <- D_boot[idx_sg, , drop = FALSE]
        DX_sg <- DX_boot[idx_sg, , drop = FALSE]

        for (j in seq_len(num_combos)) {
          if (!is.null(select) && select[i, s, j] != 1) next

          t1 <- combo[j, 1]
          t2 <- combo[j, 2]

          keep <- D_sg[,1] %in% c(t1, t2)

          if (sum(keep) < 2) next

          treat <- as.integer(D_sg[keep, 1] == t2)
          cur_X <- DX_sg[keep, , drop = FALSE]
          cur_X[,1] <- treat
          cur_y <- Y_sg[keep, i]

          # Use original subgroup but bootstrap treatment (like Stata)
          pi_2 <- sum(subgroup == sg & D_boot[,1] == t2) / n
          pi_1 <- sum(subgroup == sg & D_boot[,1] == t1) / n
          pi_z <- sum(subgroup == sg) / n

          est_b <- runreg(cur_X, cur_y, Xbar_b, Xvar_b, pi_2, pi_1, pi_z)

          # Fix test statistic to match Stata
          stat_b <- abs(est_b["ATE"] - observed_coef[i,s,j]) / (if (studentized) est_b["SE"] else 1)

          boot_stats[b,i,s,j] <- stat_b

        }
      }
    }

    if (b %% 1000 == 0) {
      cat("---------------------------------------------------------------------------------", b, "\n")
    }
  }


  return(list(stat = observed_stat,
              coef = observed_coef,
              boot = boot_stats))
}
