#' Run adjusted regression to estimate treatment effect
#'
#' @param X     A matrix: first column is treatment dummy (1 for treatment, 0 for control), rest are controls (optional)
#' @param y     A numeric vector of outcomes
#' @param Xbar  A vector of covariate means for the subgroup (can be 0)
#' @param Xvar  A covariance matrix of covariates (can be 0)
#' @param pi_2  Proportion of treatment group
#' @param pi_1  Proportion of control group
#' @param pi_z  Proportion of the whole subgroup
#'
#' @return      A named vector with ATE and SE
#' @export
runreg <- function(X, y, Xbar, Xvar, pi_2, pi_1, pi_z)  {

  D <- X[, 1]
  covariates <- if (ncol(X) > 1) X[, -1, drop = FALSE] else matrix(0, nrow=length(y), ncol=0)
  n <- nrow(covariates)

  D1_idx <- D == 1
  D0_idx <- D == 0

  y1 <- y[D1_idx]
  y0 <- y[D0_idx]

  n1 <- length(y1)
  n0 <- length(y0)

  # With covariates
  if (ncol(X) > 1) {
    X1 <- X[D1_idx, -1, drop = FALSE]
    X0 <- X[D0_idx, -1, drop = FALSE]
    DX1 <- cbind(1, X1)
    DX0 <- cbind(1, X0)
  } else {
    X1 <- matrix(0, nrow = n1, ncol = 0)
    X0 <- matrix(0, nrow = n0, ncol = 0)
    DX1 <- matrix(1, nrow = n1, ncol = 1)
    DX0 <- matrix(1, nrow = n0, ncol = 1)
  }

  # Run separate regressions
  b1 <- solve(t(DX1) %*% DX1) %*% t(DX1) %*% y1
  b0 <- solve(t(DX0) %*% DX0) %*% t(DX0) %*% y0

  if (ncol(X) > 1) {
    bX1 <- b1[-1]
    bX0 <- b0[-1]

    # Use Stata's quadvariance approach instead of var()
    e1 <- y1 - X1 %*% bX1 - b1[1]
    e0 <- y0 - X0 %*% bX0 - b0[1]

    # This matches Stata's quadvariance calculation
    s1 <- sum(e1^2) / (n1-1)
    s0 <- sum(e0^2) / (n0-1)
  } else {
    bX1 <- bX0 <- 0

    # Again using quadvariance approach
    s1 <- sum((y1 - mean(y1))^2) / (n1-1)
    s0 <- sum((y0 - mean(y0))^2) / (n0-1)
  }

  # Adjusted treatment effect estimate
  ATE <- (b1[1] - b0[1]) + sum(Xbar * (bX1 - bX0))

  # Use Stata's exact approach for the covariance term
  cov_term <- if (ncol(X) > 1) as.numeric(t(bX1 - bX0) %*% Xvar %*% (bX1 - bX0)) else 0

  # Compute SE exactly as Stata does
  SE <- sqrt((1 / pi_2) * s1 + (1 / pi_1) * s0 + (1 / pi_z) * cov_term)

  return(c(ATE = ATE, SE = SE))
}
