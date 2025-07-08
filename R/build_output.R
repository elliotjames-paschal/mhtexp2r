#' Build output data frame with multiple testing corrections
#'
#' @param stat 3D array of test statistics (outcomes x subgroups x comparisons)
#' @param coef 3D array of treatment effect coefficients (outcomes x subgroups x comparisons)
#' @param combo Matrix of treatment comparison pairs (numpc x 2)
#' @param alpha_sin 3D array of single hypothesis adjusted p-values (Remark 3.2)
#' @param alpha_mul Vector of multiple hypothesis adjusted p-values (Theorem 3.1)
#' @param pvals 3D array of p-values (1-p format)
#' @param alpha_mulm Vector of transitivity-corrected p-values (Remark 3.8)
#' @param select 3D array indicating which hypotheses to include (default: all)
#'
#' @return Data frame with hypothesis test results and corrections
#' @export
build_output <- function(stat, coef, combo, alpha_sin, alpha_mul, pvals, alpha_mulm, select = NULL) {
  dims <- dim(stat)
  num_outcomes <- dims[1]
  num_subgroups <- dims[2]
  num_comparisons <- dims[3]

  out_list <- list()
  result_counter <- 1  # Index into the result vectors

  for (i in 1:num_outcomes) {
    for (j in 1:num_subgroups) {
      for (k in 1:num_comparisons) {
        # Only include if selected (like Stata)
        if (is.null(select) || select[i, j, k] == 1) {
          out_list[[result_counter]] <- data.frame(
            outcome = i,
            subgroup = j,
            comparison = k,
            t1 = combo[k, 1],
            t2 = combo[k, 2],
            coefficient = coef[i, j, k],
            test_stat = stat[i, j, k],
            Remark3_2 = alpha_sin[i, j, k],
            Thm3_1 = alpha_mul[result_counter],
            Remark3_8 = alpha_mulm[result_counter]
          )
          result_counter <- result_counter + 1
        }
      }
    }
  }

  out <- do.call(rbind, out_list)

  pvec <- out$Remark3_2
  nh <- length(pvec)

  out$Bonf <- pmin(1, pvec * nh)

  out$Holm <- {
    o <- order(pvec)
    pvec_sorted <- pvec[o]
    holm_adjust <- pvec_sorted * rev(seq_len(nh))
    holm_adjust <- pmin(cummax(holm_adjust), 1)
    result <- numeric(nh)
    result[o] <- holm_adjust
    result
  }

  return(out)
}
