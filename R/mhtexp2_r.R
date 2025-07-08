
#' mhtexp2_r: R wrapper for multiple hypothesis testing procedure
#'
#' @param Y Matrix or data.frame of outcomes
#' @param treatment Vector of treatment assignments
#' @param controls Optional matrix or data.frame of controls
#' @param subgroup Optional vector of subgroup identifiers
#' @param combo "pairwise" or "treatmentcontrol"
#' @param bootstrap Number of bootstrap replications
#' @param studentized Logical, whether to studentize test stats (default = TRUE)
#' @param transitivity_check Logical, whether to apply transitivity correction (Remark 3.8)
#'
#' @return A list containing the output data.frame and intermediate arrays
#' @export
mhtexp2_r <- function(Y, treatment, controls = NULL, subgroup = NULL,
                      combo = "treatmentcontrol", bootstrap = 3000,
                      studentized = TRUE, transitivity_check = TRUE) {
  # Prepare data
  Y <- as.matrix(Y)
  D <- matrix(as.numeric(treatment), ncol = 1)
  X <- if (!is.null(controls)) as.matrix(controls) else NULL
  subgroup <- if (is.null(subgroup)) rep(1, nrow(Y)) else subgroup

  combo_mat <- build_combo(groups = sort(unique(D)), method = combo)

  numoc <- ncol(Y)
  numsub <- length(unique(subgroup))
  numpc <- nrow(combo_mat)

  select <- array(1, dim = c(numoc, numsub, numpc))

  # Generate random bootstrap samples
  idbootmat <- replicate(bootstrap, sample(seq_len(nrow(Y)), replace = TRUE))

  # Run bootstrap
  results <- bootstrap_runreg(
    Y = Y,
    D = D,
    DX = cbind(D, X),
    combo = combo_mat,
    B = bootstrap,
    subgroup = subgroup,
    studentized = studentized,
    idbootmat = idbootmat,
    select = select
  )

  # Build bootstrap p-value array for threshold routines
  raw_boot <- results$boot
  dims_boot <- dim(raw_boot)
  pboot <- array(NA, dim = dims_boot)
  B <- dims_boot[1]
  for (i in seq_len(dims_boot[2])) {
    for (j in seq_len(dims_boot[3])) {
      for (k in seq_len(dims_boot[4])) {
        tmp <- raw_boot[, i, j, k]
        for (l in 1:B) {
          pboot[l, i, j, k] <- 1 - (sum(tmp >= tmp[l]) / B)
        }
      }
    }
  }


  # Calculate p-values and threshold corrections
  pvals  <- calculate_pvals(results$stat, results$boot)
  alpha1 <- calculate_alphasin(pvals, pboot)

  thresholds <- calculate_alpha_unified(
    pvals = pvals,
    coefficients = results$coef,
    pboot = pboot,
    combo = combo_mat,
    alphasin = alpha1,
    select = select,
    transitivitycheck = transitivity_check
  )

  alpha2 <- thresholds$alphamul
  alpha3 <- thresholds$alphamulm

  # Return results
  return(list(
    output      = build_output(
      stat       = results$stat,
      coef       = results$coef,
      combo      = combo_mat,
      alpha_sin  = alpha1,
      alpha_mul  = alpha2,
      pvals      = pvals,
      alpha_mulm = alpha3,
      select     = select
    ),
    stat        = results$stat,
    coef        = results$coef,
    boot        = results$boot,
    pboot       = pboot,
    pvals       = pvals,
    alpha_sin   = alpha1,
    alpha_mul   = alpha2,
    alpha_mulm  = alpha3
  ))
}
