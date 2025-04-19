#' Bootstrap adjusted regression with subgroup support
#'
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
  n <- nrow(Y)
  num_outcomes <- ncol(Y)
  num_combos <- nrow(combo)

  if (is.null(subgroup)) {
    subgroup <- rep(1, n)
  }

  subgroups <- sort(unique(subgroup))
  num_subgroups <- length(subgroups)

  observed_stat <- array(NA, dim = c(num_outcomes, num_subgroups, num_combos))
  observed_coef <- array(NA, dim = c(num_outcomes, num_subgroups, num_combos))
  boot_stats     <- array(NA, dim = c(B, num_outcomes, num_subgroups, num_combos))

  for (s in seq_along(subgroups)) {
    sg <- subgroups[s]
    idx_sg <- which(subgroup == sg)

    Y_sg <- Y[idx_sg, , drop = FALSE]
    D_sg <- D[idx_sg, , drop = FALSE]
    DX_sg <- DX[idx_sg, , drop = FALSE]

    for (i in seq_len(num_outcomes)) {
      yi <- Y_sg[, i]
      for (j in seq_len(num_combos)) {

        if (!is.null(select) && select[i, s, j] != 1) next

        t1 <- combo[j, 1]
        t2 <- combo[j, 2]
        keep <- D_sg[, 1] %in% c(t1, t2)
        if (sum(keep) < 2) next

        treat <- as.integer(D_sg[keep, 1] == t2)
        cur_X <- DX_sg[keep, , drop = FALSE]
        cur_X[, 1] <- treat
        cur_y <- yi[keep]

        controls <- if (ncol(DX_sg) > 1) cur_X[, -1, drop = FALSE] else NULL
        Xbar_sg <- if (!is.null(controls)) colMeans(controls) else 0
        Xvar_sg <- if (!is.null(controls)) stats::var(controls) else 0

        est <- runreg(cur_X, cur_y, Xbar_sg, Xvar_sg,
                      mean(treat == 1), mean(treat == 0), 1)
        stat <- abs(est["ATE"]) / (if (studentized) est["SE"] else 1)

        observed_coef[i, s, j] <- est["ATE"]
        observed_stat[i, s, j] <- stat

        for (b in seq_len(B)) {
          boot_idx <- idbootmat[, b]
          idx <- boot_idx[boot_idx %in% idx_sg]
          if (length(idx) < 2) next

          yb  <- Y[idx, i]
          Db  <- D[idx, , drop = FALSE]
          DXb <- DX[idx, , drop = FALSE]

          keep_b <- Db[, 1] %in% c(t1, t2)
          if (sum(keep_b) < 2) next

          treat_b  <- as.integer(Db[keep_b, 1] == t2)
          cur_Xb   <- DXb[keep_b, , drop = FALSE]
          cur_Xb[, 1] <- treat_b
          cur_yb   <- yb[keep_b]

          controls_b <- if (ncol(DXb) > 1) cur_Xb[, -1, drop = FALSE] else NULL
          Xbar_b     <- if (!is.null(controls_b)) colMeans(controls_b) else 0
          Xvar_b     <- if (!is.null(controls_b)) stats::var(controls_b) else 0

          est_b <- runreg(cur_Xb, cur_yb, Xbar_b, Xvar_b,
                          mean(treat_b == 1), mean(treat_b == 0), 1)
          stat_b <- abs(est_b["ATE"]) / (if (studentized) est_b["SE"] else 1)
          boot_stats[b, i, s, j] <- stat_b
        }
      }
    }
  }

  return(list(
    stat = observed_stat,
    coef = observed_coef,
    boot = boot_stats
  ))
}
