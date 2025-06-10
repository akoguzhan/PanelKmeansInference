#' Grid SU p-value combination test using median
#'
#' Combines dependent p-values using the inverse moment method over a grid of r-values,
#' and aggregates the resulting pseudo p-values via median (more stable than min).
#'
#' @param pvals Numeric vector of p-values (in (0,1)).
#' @param grid_min Minimum r value in the grid (default: 5)
#' @param grid_max Maximum r value in the grid (default: 50)
#' @param n_grid Number of grid points (default: 20)
#' @param adjust Logical; if TRUE, clip p-values to (1e-10, 1 - 1e-10) to avoid instability.
#'
#' @return A list with median p-value, grid of r values, and vector of all pseudo p-values.
#' @export
grid_iu_test_median <- function(pvals, grid_min = 5, grid_max = 50, n_grid = 20, adjust = TRUE) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  if (adjust) pvals <- pmin(pmax(pvals, 1e-10), 1 - 1e-10)

  r_grid <- seq(grid_min, grid_max, length.out = n_grid)
  n <- length(pvals)

  pseudo_pvals <- sapply(r_grid, function(r) {
    SU_stat <- (mean(pvals^(-r)))^(1 / r)
    p_su <- min(r / (r - 1) / SU_stat, 1)
    return(p_su)
  })

  median_p <- median(pseudo_pvals)

  return(list(
    adjusted_p = median_p,
    r_grid = r_grid,
    pseudo_pvals = pseudo_pvals
  ))
}

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
#'
#' @return A list with:
#' \describe{
#'   \item{p_grid}{Vector of pseudo p-values across the r-grid}
#'   \item{r_grid}{The grid of r values}
#'   \item{min_p}{Minimum pseudo p-value}
#'   \item{selected_r}{The r value achieving min_p}
#' }
#' @export
grid_iu_test <- function(pvals, r_min = 5, r_max = 50, n_grid = 25) {
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

  list(
    p_grid = p_grid,
    r_grid = r_grid,
    min_p = min_p,
    selected_r = r_star
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
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
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
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }

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

#' Cauchy Combination Test (CCT)
#'
#' Combines dependent p-values using the Cauchy combination test proposed by Liu & Xie (2020).
#'
#' @param p numeric vector of p-values (all in (0,1))
#' @param weights optional numeric vector of weights (must sum to 1); default: equal weights
#'
#' @return Combined p-value
#' @export
cauchy_combine <- function(p, weights = NULL) {
  if (any(p < 0 | p > 1)) stop("All p-values must lie strictly in [0,1].")
  
  n <- length(p)
  if (is.null(weights)) {
    weights <- rep(1 / n, n)
  } else {
    if (length(weights) != n || any(weights < 0)) stop("Weights must be non-negative and match p-values in length.")
    weights <- weights / sum(weights)  # normalize
  }

  T_stat <- sum(weights * tan((0.5 - p) * pi))
  p_combined <- 0.5 - atan(T_stat) / pi
  return(min(max(p_combined, 0), 1))  # ensure in [0,1]
}

#' Combined Selective p-value Using Harmonic Mean
#'
#' Computes the combined p-value for testing the separation between two clusters
#' using the harmonic mean of selective p-values between adjacent clusters.
#' This corresponds to the formula at the bottom of page 5 in Hivert et al. (2024).
#'
#' @param pvals A numeric vector of selective p-values between adjacent clusters (length M - 1).
#'
#' @return A combined p-value between 0 and 1.
#' @references Hivert et al. (2024). Post-clustering difference testing. CSDA.
#' @export
combined_selective_harmonic <- function(pvals) {
  M <- length(pvals) + 1
  harmonic_mean <- (M - 1) / sum(1 / pvals)
  p_comb <- min(exp(1) * log(M - 1) * harmonic_mean, 1)
  return(p_comb)
}