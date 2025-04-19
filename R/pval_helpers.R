#' Calculate p-values from observed and bootstrap test statistics
#'
#' @param observed 3D array (outcomes x subgroups x comparisons)
#' @param boot 4D array (B x outcomes x subgroups x comparisons)
#'
#' @return 3D array of p-values (same shape as observed)
#' @export
calculate_pvals <- function(observed, boot) {
  B <- dim(boot)[1]
  out <- array(NA, dim = dim(observed))

  for (i in seq_len(dim(observed)[1])) {
    for (j in seq_len(dim(observed)[2])) {
      for (k in seq_len(dim(observed)[3])) {
        stat_obs <- observed[i, j, k]
        stat_boot <- boot[, i, j, k]
        out[i, j, k] <- mean(stat_boot >= stat_obs)
      }
    }
  }

  return(out)
}

#' Calculate single hypothesis adjusted thresholds (Remark 3.2)
#'
#' @param observed 3D array (outcomes x subgroups x comparisons)
#' @param boot 4D array (B x outcomes x subgroups x comparisons)
#'
#' @return 3D array of adjusted alpha thresholds (same shape)
#' @export
calculate_alphasin <- function(observed, boot) {
  B <- dim(boot)[1]
  out <- array(NA, dim = dim(observed))

  for (i in seq_len(dim(observed)[1])) {
    for (j in seq_len(dim(observed)[2])) {
      for (k in seq_len(dim(observed)[3])) {
        obs_stat <- observed[i, j, k]
        boot_stats <- boot[, i, j, k]
        sorted_boot <- sort(boot_stats, decreasing = TRUE)
        v <- obs_stat >= sorted_boot
        first_idx <- which(v)[1]
        out[i, j, k] <- if (is.na(first_idx)) 1 else first_idx / B
      }
    }
  }

  return(out)
}
