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
  covariates <- X[, -1, drop = FALSE]

  Xbar <- colMeans(covariates)
  Xvar <- stats::var(covariates)
  D1_idx <- D == 1
  D0_idx <- D == 0

  y1 <- y[D1_idx]
  y0 <- y[D0_idx]

  # With covariates
  if (ncol(X) > 1) {
    X1 <- X[D1_idx, -1, drop = FALSE]
    X0 <- X[D0_idx, -1, drop = FALSE]
    DX1 <- cbind(1, X1)
    DX0 <- cbind(1, X0)
  } else {
    X1 <- matrix(0, nrow = length(y1), ncol = 0)
    X0 <- matrix(0, nrow = length(y0), ncol = 0)
    DX1 <- matrix(1, nrow = length(y1), ncol = 1)
    DX0 <- matrix(1, nrow = length(y0), ncol = 1)
  }

  # Run separate regressions
  b1 <- solve(t(DX1) %*% DX1) %*% t(DX1) %*% y1
  b0 <- solve(t(DX0) %*% DX0) %*% t(DX0) %*% y0

  if (ncol(X) > 1) {
    bX1 <- b1[-1]
    bX0 <- b0[-1]
    s1 <- var(y1 - X1 %*% bX1)
    s0 <- var(y0 - X0 %*% bX0)
  } else {
    bX1 <- bX0 <- 0
    s1 <- var(y1)
    s0 <- var(y0)
  }

  # Adjusted treatment effect estimate
  ATE <- (b1[1] - b0[1]) + sum(Xbar * (bX1 - bX0))

  # Adjusted SE
  SE <- sqrt((1 / pi_2) * s1 + (1 / pi_1) * s0 + (1 / pi_z) * sum((bX1 - bX0)^2 * diag(Xvar)))

  return(c(ATE = ATE, SE = SE))
}
