# Main user-facing function
# File: mhtexp2_r.R

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
                      studentized = TRUE, transitivity_check = TRUE,
                      exclude = NULL, only = NULL,
                      idbootmat = NULL, treatnames = NULL) {
  # Prepare data
  Y <- as.matrix(Y)
  D <- matrix(as.numeric(treatment), ncol = 1)
  X <- if (!is.null(controls)) as.matrix(controls) else NULL
  subgroup <- if (is.null(subgroup)) rep(1, nrow(Y)) else subgroup

  # Build combo, select, and idbootmat
  combo_mat <- build_combo(groups = sort(unique(D)), method = combo)
  numoc <- ncol(Y)
  numsub <- length(unique(subgroup))
  numpc <- nrow(combo_mat)

  select <- array(1, dim = c(numoc, numsub, numpc))
  if (!is.null(only)) {
    select[] <- 0
    for (row in seq_len(nrow(only))) {
      select[only[row,1], only[row,2], only[row,3]] <- 1
    }
  }
  if (!is.null(exclude)) {
    for (row in seq_len(nrow(exclude))) {
      select[exclude[row,1], exclude[row,2], exclude[row,3]] <- 0
    }
  }

  if (is.null(idbootmat)) {
    idbootmat <- replicate(bootstrap, sample(seq_len(nrow(Y)), replace = TRUE))
  }

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

  # Build bootstrap p-value array (pboot) for threshold routines
  raw_boot <- results$boot
  dims_boot <- dim(raw_boot)
  pboot <- array(NA, dim = dims_boot)
  B <- dims_boot[1]
  for (i in seq_len(dims_boot[2])) {
    for (j in seq_len(dims_boot[3])) {
      for (k in seq_len(dims_boot[4])) {
        tmp <- raw_boot[, i, j, k]
        pboot[, i, j, k] <- sapply(tmp, function(x) 1 - mean(tmp >= x))
      }
    }
  }

  # P-value + threshold corrections
  pvals  <- calculate_pvals(results$stat, results$boot)
  alpha1 <- calculate_alphasin(pvals, pboot)
  alpha2 <- calculate_alphamul(pvals, pboot, alpha1)

  # Build ID vectors
  dim_obs <- dim(results$stat)
  nh <- prod(dim_obs)
  outcome_ids  <- rep(seq_len(dim_obs[1]), each = dim_obs[2] * dim_obs[3])
  subgroup_ids <- rep(rep(seq_len(dim_obs[2]), each = dim_obs[3]), times = dim_obs[1])
  combo_ids    <- rep(seq_len(dim_obs[3]), times = dim_obs[1] * dim_obs[2])
  combo_long   <- combo_mat[combo_ids, , drop = FALSE]

  # Optional transitivity correction (Remark 3.8)
  alpha3 <- if (transitivity_check) {
    calculate_alphamulm(
      observed = pvals,     # This is already correct!
      boot     = pboot,
      combo    = combo_long,
      outcome_ids  = outcome_ids,
      subgroup_ids = subgroup_ids
    )
  } else {
    alpha2
  }

  # Return results
  return(list(
    output      = build_output(
      stat       = results$stat,
      coef       = results$coef,
      combo      = combo_mat,
      alpha_sin  = alpha1,
      alpha_mul  = alpha2,
      pvals      = pvals,
      alpha_mulm = alpha3
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


