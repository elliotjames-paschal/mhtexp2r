#' Calculate multiple-testing adjusted thresholds (Theorem 3.1)
#' Calculate multiple-testing adjusted thresholds (Theorem 3.1) with subgroup support
#'
#' @param observed 3D array of observed test stats (outcome x subgroup x comparison)
#' @param boot 4D array of bootstrap stats (B x outcome x subgroup x comparison)
#'
#' @return Matrix of adjusted alpha thresholds (same shape as observed)
#' @export
calculate_alphamul <- function(observed, boot) {
  B <- dim(boot)[1]
  dim_obs <- dim(observed)
  nh <- prod(dim_obs)

  obs_vec <- as.vector(observed)
  boot_mat <- matrix(NA, nrow = B, ncol = nh)

  idx <- 1
  for (i in seq_len(dim_obs[1])) {
    for (j in seq_len(dim_obs[2])) {
      for (k in seq_len(dim_obs[3])) {
        boot_mat[, idx] <- boot[, i, j, k]
        idx <- idx + 1
      }
    }
  }

  ord <- order(obs_vec)
  obs_vec <- obs_vec[ord]
  boot_mat <- boot_mat[, ord, drop = FALSE]

  alphamul <- numeric(nh)
  for (i in seq_len(nh)) {
    maxstats <- apply(boot_mat[, i:nh, drop = FALSE], 1, max)
    sort_max <- sort(maxstats, decreasing = TRUE)
    exceed <- which(obs_vec[i] >= sort_max)
    alphamul[i] <- if (length(exceed) == 0) 1 else min(exceed) / B
  }

  alphamul_out <- numeric(nh)
  alphamul_out[ord] <- alphamul
  array(alphamul_out, dim = dim_obs)
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
calculate_alphamulm <- function(observed, boot, combo, outcome_ids, subgroup_ids) {
  B <- dim(boot)[1]
  nh <- length(outcome_ids)

  obs_vec <- as.vector(observed)
  boot_mat <- matrix(NA, nrow = B, ncol = nh)

  idx <- 1
  for (i in seq_len(dim(observed)[1])) {
    for (j in seq_len(dim(observed)[2])) {
      for (k in seq_len(dim(observed)[3])) {
        boot_mat[, idx] <- boot[, i, j, k]
        idx <- idx + 1
      }
    }
  }

  if (length(obs_vec) != ncol(boot_mat) ||
      length(obs_vec) != nrow(combo) ||
      length(obs_vec) != length(outcome_ids) ||
      length(obs_vec) != length(subgroup_ids)) {
    stop("Mismatch in dimensions: check combo, outcome_ids, subgroup_ids.")
  }

  ord <- order(obs_vec)
  obs_vec <- obs_vec[ord]
  boot_mat <- boot_mat[, ord, drop = FALSE]
  combo <- combo[ord, , drop = FALSE]
  outcome_ids <- outcome_ids[ord]
  subgroup_ids <- subgroup_ids[ord]

  alphamulm <- numeric(nh)
  rejected <- list()

  for (i in seq_len(nh)) {
    if (i == 1) {
      alphamulm[i] <- 1
      rejected[[1]] <- list(outcome = outcome_ids[i], subgroup = subgroup_ids[i],
                            t1 = combo[i, 1], t2 = combo[i, 2])
      next
    }

    found <- FALSE
    for (j in seq(nh - i + 1, 1)) {
      subset_indices <- combn(i:nh, j, simplify = FALSE)
      for (subset in subset_indices) {
        subset_hypotheses <- lapply(subset, function(k) {
          list(outcome = outcome_ids[k], subgroup = subgroup_ids[k],
               t1 = combo[k, 1], t2 = combo[k, 2])
        })

        if (!violates_transitivity(subset_hypotheses, rejected)) {
          maxstats <- apply(boot_mat[, unlist(subset), drop = FALSE], 1, max)
          sort_max <- sort(maxstats, decreasing = TRUE)
          q <- which(obs_vec[i] >= sort_max)
          alphamulm[i] <- if (length(q) == 0) 1 else min(q) / B

          rejected[[length(rejected) + 1]] <- list(outcome = outcome_ids[i],
                                                   subgroup = subgroup_ids[i],
                                                   t1 = combo[i, 1],
                                                   t2 = combo[i, 2])
          found <- TRUE
          break
        }
      }
      if (found) break
    }

    if (!found) alphamulm[i] <- 1
  }

  alphamulm_out <- numeric(nh)
  alphamulm_out[ord] <- alphamulm
  array(alphamulm_out, dim = dim(observed))
}
