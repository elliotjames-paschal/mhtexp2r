{

  dims <- dim(stat)
  num_outcomes <- dims[1]
  num_subgroups <- dims[2]
  num_comparisons <- dims[3]

  # Preserve original row order
  out <- expand.grid(
    outcome = seq_len(num_outcomes),
    subgroup = seq_len(num_subgroups),
    comparison = seq_len(num_comparisons)
  )

  # Add results
  out$coefficient <- as.vector(coef)
  out$test_stat   <- as.vector(stat)
  out$Remark3_2   <- as.vector(alpha_sin)
  out$Thm3_1      <- as.vector(alpha_mul)
  out$Remark3_8   <- as.vector(alpha_mulm)

  # Add Bonferroni and Holm (based on single-test p-values)
  pvec <- out$Remark3_2
  nh <- length(pvec)

  # Bonferroni correction
  out$Bonf <- pmin(1, pvec * nh)
  # Bonferroni correction
  bonf_flat <- pmin(1, alphasin_flat * num_hypotheses)

  # Holm correction
  out$Holm <- {
    o <- order(pvec)
    pvec_sorted <- pvec[o]

    # Holm adjustment
    holm_adjust <- pvec_sorted * rev(seq_len(nh))
    holm_adjust <- pmin(cummax(holm_adjust), 1)

    result <- numeric(nh)
    result[o] <- holm_adjust
    result
  }

  # Add treatment-control pairs
  combo_rep <- combo[rep(seq_len(nrow(combo)), each = num_outcomes * num_subgroups), ]
  out$t1 <- combo_rep[, 1]
  out$t2 <- combo_rep[, 2]

  # Final formatting
  out$comparison <- NULL  # Remove unneeded column

  out <- out[, c("outcome", "subgroup", "t1", "t2",
                 "coefficient", "Remark3_2", "Thm3_1", "Remark3_8",
                 "Bonf", "Holm")]

  return(out)
}
