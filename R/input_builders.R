#' Build Y matrix (outcomes)
#' @param data Data frame
#' @param outcomes Character vector of outcome column names
#' @return A matrix of outcomes
#' @export
build_Y <- function(data, outcomes) {
  as.matrix(data[, outcomes, drop = FALSE])
}

#' Build D matrix (treatment indicator)
#' @param data Data frame
#' @param treatment Character vector of outcome column names
#' @return A matrix of outcomes
#' @export
build_D <- function(data, treatment) {
  as.matrix(data[, treatment, drop = FALSE])
}

#' Build DX matrix (treatment + controls)
#' @param data Data frame
#' @param treatment Character vector of outcome column names
#' @param controls Optional character vector of control variable names
#' @return A matrix of outcomes
#' @export
build_DX <- function(data, treatment, controls = NULL) {
  vars <- c(treatment, controls)
  as.matrix(data[, vars, drop = FALSE])
}

#' Build subgroup vector
#' @param data Data frame
#' @param subgroup Optional character name of subgroup column
#' @param n Integer, number of rows
#' @return Integer vector of subgroup IDs
#' @export
build_sub <- function(data, subgroup, n = nrow(data)) {
  if (is.null(subgroup) || subgroup == "") {
    rep(1, n)
  } else {
    as.integer(data[[subgroup]])
  }
}

#' Build treatment-control combinations
#' @param groups A vector of unique treatment groups (assumed sorted and numeric)
#' @param method String: either "pairwise" or "treatmentcontrol"
#'
#' @return A matrix of treatment-control combinations
#' @export
build_combo <- function(groups, method = "treatmentcontrol") {
  groups <- sort(unique(groups))

  if (method == "pairwise") {
    combo <- t(combn(groups, 2))
  } else if (method == "treatmentcontrol") {
    control <- min(groups)
    combo <- cbind(rep(control, length(groups) - 1), groups[groups != control])
  } else {
    stop("method must be 'pairwise' or 'treatmentcontrol'")
  }

  return(combo)
}
