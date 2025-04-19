#' Format final output from mhtexp2_r
#'
#' @param observed 3D array of observed stats (outcome x subgroup x comparison)
#' @param combo Matrix of treatment-control pairs
#' @param alpha_sin Matrix of single-hypothesis thresholds
#' @param alpha_mul Matrix of multiple-testing corrected thresholds
#' @param pvals Matrix of p-values
#' @param alpha_mulm Matrix of transitivity-corrected thresholds
#'
#' @return A long-format data.frame of test results
#' @export
build_output <- function(observed, combo, alpha_sin, alpha_mul, pvals, alpha_mulm) {
  dims <- dim(observed)
  num_outcomes <- dims[1]
  num_subgroups <- dims[2]
  num_comparisons <- dims[3]

  # Create flat index
  out <- expand.grid(
    outcome = seq_len(num_outcomes),
    subgroup = seq_len(num_subgroups),
    comparison = seq_len(num_comparisons)
  )

  # Add stats
  out$stat <- as.vector(observed)
  out$pval <- as.vector(pvals)
  out$alpha_sin <- as.vector(alpha_sin)
  out$alpha_mul <- as.vector(alpha_mul)
  out$alpha_mulm <- as.vector(alpha_mulm)

  # Add treatment-control pairs
  combo_rep <- combo[rep(seq_len(nrow(combo)), each = num_outcomes * num_subgroups), ]
  out$t1 <- combo_rep[, 1]
  out$t2 <- combo_rep[, 2]

  return(out)
}
