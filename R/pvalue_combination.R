
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