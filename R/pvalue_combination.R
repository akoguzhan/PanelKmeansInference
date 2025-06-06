#' Grid Inverse-U Test (Adaptive Spreng–Urga Combination)
#'
#' Combines multiple p-values using an adaptive inverse-moment (Inverse-U) rule
#' proposed by Spreng and Urga (2023), evaluated over a grid of inverse powers.
#' This method is powerful under sparse alternatives and robustified via grid search.
#'
#' @param pvals A numeric vector of p-values (all in [0, 1]).
#' @param r_min Minimum inverse moment power (default: 5).
#' @param r_max Maximum inverse moment power (default: 50).
#' @param n_grid Number of grid points (default: 25).
#' @param adjust Logical; whether to apply Bonferroni correction for grid search (default: TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{p_grid}{Vector of pseudo p-values across the r-grid}
#'   \item{r_grid}{The grid of r values}
#'   \item{min_p}{Minimum pseudo p-value}
#'   \item{selected_r}{The r value achieving min_p}
#'   \item{adjusted_p}{Bonferroni-adjusted combined p-value if \code{adjust = TRUE}}
#' }
#' @export
grid_iu_test <- function(pvals, r_min = 5, r_max = 50, n_grid = 25, adjust = TRUE) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n <- length(pvals)
  r_grid <- seq(r_min, r_max, length.out = n_grid)

  # Vectorized calculation
  mean_pvals_neg_r <- vapply(r_grid, function(r) mean(pvals^(-r)), numeric(1))
  F_r <- mean_pvals_neg_r^(1 / r_grid)
  p_grid <- pmin((r_grid / (r_grid - 1)) * (1 / F_r), 1)

  min_idx <- which.min(p_grid)
  min_p <- p_grid[min_idx]
  r_star <- r_grid[min_idx]
  adj_p <- if (adjust) min(1, n_grid * min_p) else min_p

  list(
    p_grid = p_grid,
    r_grid = r_grid,
    min_p = min_p,
    selected_r = r_star,
    adjusted_p = adj_p
  )
}

#' IU P-value Combination Test (Spreng and Urga, 2022)
#'
#' Combines a vector of p-values using the Inverse-Uniform (IU) method.
#'
#' @param pvals A numeric vector of p-values (all in (0, 1]).
#' @param r Power parameter (default = 20).
#' @param alpha Significance level (default = 0.05).
#'
#' @return A list with test statistic, critical value, rejection indicator, and combined p-value.
#' @references Spreng, L., & Urga, G. (2022). Combining p-values for Multivariate Predictive Ability Testing.
#' @export
iu_test <- function(pvals, r = 20, alpha = 0.05) {
  stopifnot(is.numeric(pvals), length(pvals) > 0, all(pvals > 0 & pvals <= 1))
  n <- length(pvals)
  Prn <- (mean(pvals^(-r)))^(1 / r)
  critical_value <- r / (alpha * (r - 1))
  reject <- Prn > critical_value
  combined_p <- min(r / ((r - 1) * Prn), 1)
  list(
    statistic = Prn,
    critical_value = critical_value,
    reject = reject,
    pval_combined = combined_p
  )
}

#' Grid Harmonic P-Value Combination (Vovk et al., 2022)
#'
#' This function combines a vector of p-values using the Grid Harmonic method 
#' described in Vovk, Wang, and Wang (2022, Annals of Statistics, Proposition 5).
#' It calculates multiple M-mean p-values over a geometric grid of powers and returns 
#' the minimum adjusted p-value.
#'
#' @param pvals Numeric vector of p-values in (0, 1].
#' @param r_min Minimum value of the power parameter r (default: 1.05).
#' @param r_max Maximum value of the power parameter r (default: 10).
#' @param B Number of grid points to evaluate between r_min and r_max (default: 20).
#'
#' @return A list with:
#' \describe{
#'   \item{p_combined}{Combined p-value from the grid harmonic rule.}
#'   \item{r_grid}{Grid of r values used.}
#'   \item{M_values}{M-mean values for each r.}
#'   \item{adjusted_values}{Adjusted p-values at each r (min taken).}
#' }
#' 
#' @references Vovk, V., Wang, R., & Wang, R. (2022). Admissible ways of merging
#' p-values under arbitrary dependence.
#'
#' @export
grid_harmonic_pcombine <- function(pvals, r_min = 1.05, r_max = 10, B = 20) {
  if (any(pvals <= 0 | pvals > 1)) stop("All p-values must be in (0, 1].")

  n <- length(pvals)
  r_grid <- exp(seq(log(r_min), log(r_max), length.out = B))

  M_values <- sapply(r_grid, function(r) {
    (mean(pvals^r))^(1/r)
  })

  adjusted <- (r_grid + 1)^(1 / r_grid) * M_values

  p_combined <- min(adjusted, 1)

  return(list(
    pval_combined = p_combined,
    r_grid = r_grid,
    M_values = M_values,
    adjusted_values = adjusted
  ))
}

#' Tippett’s Method for P-value Combination
#'
#' Combines multiple p-values by taking their minimum and adjusting
#' for the number of tests using Tippett’s rule.
#'
#' @param pvals A numeric vector of p-values (all in [0, 1]).
#' @return A single combined p-value.
#' @examples
#' pvals <- c(0.01, 0.20, 0.05)
#' tippett_p(pvals)
#' @export
tippett_p <- function(pvals) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n <- length(pvals)
  p_min <- min(pvals)
  combined <- 1 - (1 - p_min)^n
  return(min(combined, 1))
}

#' Bonferroni Combination Method
#'
#' Combines multiple p-values using the Bonferroni correction:
#' \( p_{\text{combined}} = n \cdot \min(p_i) \), clipped at 1.
#'
#' @param pvals A numeric vector of p-values (each in [0, 1]).
#' @return A single Bonferroni-adjusted p-value.
#' @examples
#' pvals <- c(0.01, 0.20, 0.05)
#' bonferroni_p(pvals)
#' @export
bonferroni_p <- function(pvals) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n <- length(pvals)
  combined <- n * min(pvals)
  return(min(combined, 1))
}