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

        # Get observed test statistic for this hypothesis
        observed_stat <- observed[i, j, k]

        # Get all bootstrap statistics for this hypothesis
        bootstrap_stats <- boot[, i, j, k]

        # Count how many bootstrap stats >= observed stat
        exceed_count <- sum(bootstrap_stats >= observed_stat)

        # Calculate 1 - p-value (matching Stata exactly)
        one_minus_p <- 1 - (exceed_count / B)

        out[i, j, k] <- one_minus_p
      }
    }
  }

  return(out)
}


#' Calculate single hypothesis adjusted thresholds (Remark 3.2)
#'
#' @param pact 3D array (outcomes x subgroups x comparisons) - observed "1-p" values
#' @param pboot 4D array (B x outcomes x subgroups x comparisons) - bootstrap "1-p" values
#'
#' @return 3D array of adjusted alpha thresholds (same shape)
#' @export
calculate_alphasin <- function(pact, pboot) {

  B <- dim(pboot)[1]
  out <- array(NA, dim = dim(pact))

  for (i in seq_len(dim(pact)[1])) {
    for (j in seq_len(dim(pact)[2])) {
      for (k in seq_len(dim(pact)[3])) {

        # Get observed "1-p" value for this hypothesis
        observed_1p <- pact[i, j, k]

        # Get bootstrap "1-p" values for this hypothesis
        boot_1p <- pboot[, i, j, k]

        # Sort bootstrap values in descending order (largest first)
        sorted_boot <- sort(boot_1p, decreasing = TRUE)

        # Find where observed value is >= sorted bootstrap values
        v <- observed_1p >= sorted_boot

        # Find first TRUE index (where observed >= bootstrap)
        first_idx <- which(v)[1]

        # Calculate proportion
        if (is.na(first_idx)) {
          q <- 1
        } else {
          q <- first_idx / B
        }

        out[i, j, k] <- q
      }
    }
  }

  return(out)
}
