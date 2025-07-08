#' Calculate multiple hypothesis testing thresholds with optional transitivity
#'
#' @param pvals 3D array of observed p-values (outcome x subgroup x comparison)
#' @param coefficients 3D array of treatment effect coefficients (outcome x subgroup x comparison)
#' @param pboot 4D array of bootstrap p-values (B x outcome x subgroup x comparison)
#' @param combo Matrix of treatment-control pairs
#' @param alphasin 3D array of single-hypothesis adjusted p-values (outcome x subgroup x comparison)
#' @param select 3D array indicating which hypotheses to include (optional, default = all)
#' @param transitivitycheck Logical, whether to apply transitivity correction (Remark 3.8, default = FALSE)
#'
#' @return Vector of adjusted alpha thresholds in original hypothesis order
#' @export
calculate_alpha_unified <- function(pvals, coefficients, pboot, combo, alphasin, select = NULL, transitivitycheck = FALSE) {

  # Helper functions
  find_first <- function(v) {
    idx <- which(v)
    if (length(idx) == 0) return(NULL)
    return(idx[1])
  }

  nchoosek <- function(V, K) {
    if (K <= 0) return(matrix(numeric(0), 0, max(K, 0)))
    if (length(V) < K) return(matrix(numeric(0), 0, K))
    if (K > length(V)) return(matrix(numeric(0), 0, K))

    combos <- combn(V, K)
    return(t(combos))
  }

  # Stata-equivalent asarray implementation
  create_asarray <- function() {
    list(data = list(), max_key = 0)
  }

  put_asarray <- function(arr, key, value) {
    arr$data[[as.character(key)]] <- value
    if (key > arr$max_key) arr$max_key <- key
    return(arr)
  }

  get_asarray <- function(arr, key) {
    return(arr$data[[as.character(key)]])
  }

  get_max_key <- function(arr) {
    return(arr$max_key)
  }

  # Stata's ismember function - exact replication
  ismember_stata <- function(A, B, r = 1) {
    if (r == 1) {
      # Row-wise comparison
      res <- numeric(nrow(A))
      for (i in 1:nrow(A)) {
        if (all(A[i, ] == B[i, ])) {
          res[i] <- 1
        } else {
          res[i] <- 0
        }
      }
    } else {
      # Element-wise comparison for vectors
      res <- numeric(length(A))
      for (i in 1:length(A)) {
        if (any(A[i] == B)) {
          res[i] <- 1
        } else {
          res[i] <- 0
        }
      }
    }
    return(res)
  }


  # Build statsall matrix internally
  dim_obs <- dim(pvals)
  numoc <- dim_obs[1]
  numsub <- dim_obs[2]
  numpc <- dim_obs[3]
  B <- dim(pboot)[1]

  if (is.null(select)) {
    select <- array(1, dim = dim_obs)
  }

  nh <- sum(select)
  combo_for_k <- unique(combo)
  num_treatments <- ncol(combo_for_k)
  total_cols <- 3 + num_treatments + 3 + B

  # Build statsall matrix
  statsall <- matrix(0, nrow = nh, ncol = total_cols)

  counter <- 1
  for (i in 1:numoc) {
    for (j in 1:numsub) {
      for (k in 1:numpc) {
        if (select[i, j, k] == 1) {
          statsall[counter, 1] <- counter  # Sequential hypothesis ID
          statsall[counter, 2] <- i        # outcome
          statsall[counter, 3] <- j        # subgroup
          statsall[counter, 4:(3+num_treatments)] <- combo_for_k[k, ]  # treatments
          statsall[counter, 4+num_treatments] <- coefficients[i, j, k]  # coefficient
          statsall[counter, 5+num_treatments] <- alphasin[i, j, k]      # single p-value (psin)
          statsall[counter, 6+num_treatments] <- pvals[i, j, k]         # observed p-value (pact)

          # Bootstrap p-values
          for (b in 1:B) {
            statsall[counter, 6+num_treatments+b] <- pboot[b, i, j, k]
          }
          counter <- counter + 1
        }
      }
    }
  }

  # Now do the threshold calculations
  psin_col <- 5 + num_treatments
  pact_col <- 6 + num_treatments
  bootstrap_start <- 7 + num_treatments

  # Sort by single hypothesis p-values
  ord <- order(statsall[, psin_col])
  statsrank <- statsall[ord, ]

  alphamul <- numeric(nh)   # Standard results (Theorem 3.1)
  alphamulm <- numeric(nh)  # Transitivity aware results (Theorem 3.8)

  for (i in 1:nh) {
    # ALWAYS calculate standard alphamul (Theorem 3.1)
    remaining_boot <- statsrank[i:nh, bootstrap_start:(bootstrap_start+B-1), drop = FALSE]
    maxstats <- apply(remaining_boot, 2, max)
    sortmaxstats <- sort(maxstats, decreasing = TRUE)
    v <- statsrank[i, pact_col] >= sortmaxstats
    indx <- which(v)[1]

    alphamul[i] <- ifelse(is.na(indx), 1, indx / B)

    # Calculate alphamulm based on conditions
    if (i == 1 || transitivitycheck == FALSE) {
      # Copy standard result (matches Stata: alphamulm[i] = alphamul[i])
      alphamulm[i] <- alphamul[i]
    } else {
      sortmaxstatsm <- rep(0, B)

      # Loop through subset sizes from largest to smallest
      for (j in (nh - i + 1):1) {

        # Get remaining hypothesis IDs
        remaining_indices <- statsrank[i:nh, 1]

        # Generate all subsets of size j
        if (j > 0 && length(remaining_indices) >= j) {
          subsets <- nchoosek(remaining_indices, j)
        } else {
          next
        }

        if (nrow(subsets) == 0) next

        sumcont <- 0  # Count of contradictory subsets

        # Check each subset
        for (k in 1:nrow(subsets)) {
          subset <- subsets[k, ]
          cont <- 0  # Contradiction flag

          # Find subset indices in statsall - moved outside of below function bc wtf
          subset_indices <- match(subset, statsall[, 1])

          # Check against all previously rejected hypotheses
          for (l in 1:(i-1)) {

            # Get tempA and tempB exactly like Stata
            tempA <- statsall[subset_indices, 2:3, drop = FALSE]  # outcome, subgroup
            tempB <- matrix(rep(statsrank[l, 2:3], length(subset)),
                            nrow = length(subset), byrow = TRUE)

            # Find hypotheses in subset with same outcome/subgroup as rejected hypothesis l
            ismember_result <- ismember_stata(tempA, tempB, r=1)
            sameocsub <- subset[which(ismember_result == 1)]


            # Setup transitivity groups (should be BEFORE the if-else block)
            if (length(sameocsub) >= 1) {
              sameocsub_indices <- match(sameocsub, statsall[, 1])
              treatment_pairs <- statsrank[sameocsub_indices, 4:(3+num_treatments), drop = FALSE]

              tran <- create_asarray()
              for (tp in 1:nrow(treatment_pairs)) {
                tran <- put_asarray(tran, tp, treatment_pairs[tp, ])
              }
              trantemp <- tran
            }

            if (length(sameocsub) <= 1) {
              # No transitivity constraint possible with ≤1 hypothesis
              cont <- 0

              # Process the ENTIRE subset k (like Stata does)
              subset_indices <- match(subset, statsall[, 1])
              subset_boot <- statsall[subset_indices, bootstrap_start:(bootstrap_start+B-1), drop = FALSE]

              if (length(subset) == 1) {
                maxstatsm <- as.numeric(subset_boot)
              } else {
                maxstatsm <- apply(subset_boot, 2, max)
              }
              maxstatsm_sorted <- sort(maxstatsm, decreasing = TRUE)
              sortmaxstatsm <- pmax(sortmaxstatsm, maxstatsm_sorted)

              break

            } else {
              # TRANSITIVITY ANALYSIS - Multiple matching hypotheses

              # Iterative transitive closure - exact Stata while condition
              counter_tc <- 1
              while (get_max_key(tran) > get_max_key(trantemp) || counter_tc == 1) {
                tran <- trantemp
                trantemp <- create_asarray()

                # Initialize with first group
                if (get_max_key(tran) >= 1) {
                  trantemp <- put_asarray(trantemp, 1, get_asarray(tran, 1))
                }
                counter_tc <- counter_tc + 1

                # Try to merge remaining groups (exact Stata logic)
                if (get_max_key(tran) >= 2) {
                  for (m in 2:get_max_key(tran)) {
                    belong <- 0  # Track number of groups this can merge with

                    for (n in 1:get_max_key(trantemp)) {
                      trantempn <- get_asarray(trantemp, n)
                      tranm <- get_asarray(tran, m)

                      # Exact Stata logic: uniqrows( (trantempn, tranm)' )'
                      combined_matrix <- rbind(trantempn, tranm)  # Stack vertically
                      unq <- unique(as.vector(combined_matrix))   # Get unique elements
                      total_elements <- length(trantempn) + length(tranm)

                      # Stata test: unq :< cols(trantempn) + cols(tranm)
                      if (length(unq) < total_elements) {
                        # Merge groups - store unique elements
                        trantemp <- put_asarray(trantemp, n, unq)
                        belong <- belong + 1

                        # Stata's special case handling
                        if (n == get_max_key(trantemp) && belong == 0) {
                          new_key <- get_max_key(trantemp) + 1
                          trantemp <- put_asarray(trantemp, new_key, tranm)
                        }
                      }
                    }

                    if (belong == 0) {
                      new_key <- get_max_key(trantemp) + 1
                      trantemp <- put_asarray(trantemp, new_key, get_asarray(tran, m))
                    }
                  }
                }
              }

              # Check for contradiction with rejected hypothesis l
              rejected_treatments <- statsrank[l, 4:(3+num_treatments)]

              for (p in 1:get_max_key(tran)) {
                tran_group <- get_asarray(tran, p)

                # Check if both treatments from rejected hypothesis are in same group
                if (sum(ismember_stata(rejected_treatments, tran_group, r = 0)) == 2) {
                  # Both treatments from rejected hypothesis are in same group - contradiction!
                  cont <- 1
                  break
                }
              }
            }

            if (cont == 1) break  # Contradiction found, exit loop over rejected hypotheses
          }

          sumcont <- sumcont + cont

          if (cont == 0) {

            # Valid subset - process immediately (matches Stata exactly)
            subset_indices <- match(subset, statsall[, 1])
            subset_boot <- statsall[subset_indices, bootstrap_start:(bootstrap_start+B-1), drop = FALSE]

            if (length(subset) == 1) {
              maxstats_subset <- as.numeric(subset_boot)
            } else {
              maxstats_subset <- apply(subset_boot, 2, max)
            }

            maxstats_sorted <- sort(maxstats_subset, decreasing = TRUE)
            sortmaxstatsm <- pmax(sortmaxstatsm, maxstats_sorted)

          }
        }

        # Early termination if all subsets of size j are valid
        if (sumcont == 0) {
          break  # Matches Stata's break condition exactly
        } else {
        }
      }

      # Final p-value calculation
      observed_pval <- statsrank[i, pact_col]
      v <- observed_pval >= sortmaxstatsm
      indx <- find_first(v)

      if (is.null(indx)) {
        qm <- 1
      } else {
        qm <- indx / B
      }

      alphamulm[i] <- qm
    }
  }

  result_alphamul <- numeric(nh)
  result_alphamulm <- numeric(nh)
  result_alphamul[ord] <- alphamul
  result_alphamulm[ord] <- alphamulm

  # Return both results
  return(list(
    alphamul = result_alphamul,   # Theorem 3.1
    alphamulm = result_alphamulm  # Remark 3.8
  ))
}
