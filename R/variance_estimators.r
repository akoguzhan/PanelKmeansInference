#' Long-run Variance Estimation via Eigenvalue Weighted Covariance (EWC)
#'
#' @description Estimates the long-run variance using the Eigenvalue Weighted Covariance estimator.
#'
#' @param X A TxK numeric matrix.
#' @param lrv_par An integer; number of base estimators (e.g., harmonic terms).
#' If NULL, it is selected using the CPE-optimal rule from Sun (2013).
#'
#' @return A list with the KxK estimated long-run variance matrix and the lrv_par used.
#' @keywords internal
#'
EWC <- function(X, lrv_par = NULL) {
  X <- as.matrix(X)
  Tobs <- NROW(X)
  P <- NCOL(X)

  # Center the data
  X <- scale(X, scale = FALSE)

  # If lrv_par is NULL, select it using Sun's (2013) CPE-optimal rule
  if (is.null(lrv_par)) {
    lrv_par <- min(P*Tobs^(2/3),Tobs)
    # lrv_par <- Tobs
  }

  if (lrv_par < 1 || lrv_par > Tobs) stop("lrv_par must be between 1 and Tobs.")

  # Vectorized computation of all cosine basis vectors
  j_vec <- 1:lrv_par
  t_vec <- 1:Tobs
  # Create a Tobs x lrv_par matrix of cosines
  Cmat <- outer(t_vec - 0.5, j_vec, function(t, j) cos(pi * j * t / Tobs))
  # Each column of Cmat is a cosine basis vector for a given j

  # Compute all Lj at once: (P x Tobs) %*% (Tobs x lrv_par) = (P x lrv_par)
  Lj_mat <- sqrt(2 / Tobs) * t(X) %*% Cmat  # (P x lrv_par)

  # Compute OmegaHat as the average of Lj %*% t(Lj) over all j
  OmegaHat <- Lj_mat %*% t(Lj_mat) / lrv_par

  return(list(S = OmegaHat, lrv_par = lrv_par))
}

#' Long-run Variance Estimation via Newey-West
#'
#' @description Estimates the long-run variance of a time series using the Newey-West estimator.
#'
#' @param X A TxK numeric matrix.
#' @param lrv_par An integer specifying the maximum lag order. If NULL, uses 0 (sample covariance).
#'
#' @return A KxK estimated long-run variance matrix.
#' @keywords internal
#'
NeweyWest <- function(X, lrv_par = NULL) {
  X <- as.matrix(X)
  Tobs <- NROW(X)
  P <- NCOL(X)

  if (Tobs < 2) stop("Not enough observations.")
  if (is.null(lrv_par)) lrv_par <- 0
  if (lrv_par < 0 || lrv_par >= Tobs) stop("lrv_par must be between 0 and Tobs - 1.")

  # Center the data
  X <- scale(X, scale = FALSE)

  # Sample covariance
  OmegaHat <- t(X) %*% X / Tobs

  if (lrv_par > 0) {
    for (h in 1:lrv_par) {
      # Fast lagging: X[1:(Tobs-h), ] and X[(h+1):Tobs, ]
      X1 <- X[(h + 1):Tobs, , drop = FALSE]
      X2 <- X[1:(Tobs - h), , drop = FALSE]
      gamma_h <- (t(X1) %*% X2 + t(X2) %*% X1) / Tobs
      weight <- 1 - h / (lrv_par + 1)
      OmegaHat <- OmegaHat + weight * gamma_h
    }
  }
  return(OmegaHat)
}
