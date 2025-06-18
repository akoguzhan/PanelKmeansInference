#' Hybrid Bonferroni Combination with Custom P-Value Merger
#'
#' Combines a set of p-values using both Bonferroni's method and a user-supplied
#' p-value combination function, and returns twice the minimum of the two (capped at 1).
#' This follows the hybrid strategy proposed in Vovk, Wang, and Wang (2022).
#'
#' @param pvals Numeric vector of p-values in (0, 1].
#' @param combine_fun A function that takes `pvals` and returns a combined p-value.
#'
#' @return A numeric scalar: the combined p-value using the hybrid Bonferroni rule.
#'
#' @examples
#' hybrid_bonferroni_combination(runif(5), pGmean_neg)
#'
#' @export
bonferroni_compound_pcombine <- function(pvals, combine_fun1, combine_fun2) {
  # Ensure all p-values are in (0,1]
  pvals <- pmin(pmax(pvals, .Machine$double.eps), 1)

  # Compute combined p-values
  p_comb1 <- combine_fun1(pvals)
  p_comb2 <- combine_fun2(pvals)

  # Apply hybrid rule: Bonferroni over two combinations
  pval <- min(2 * min(p_comb1, p_comb2), 1)
  pval
}

#' Geometric Mean Combination (×e)
#'
#' Combines p-values using the geometric mean and multiplies by Euler's number (e ≈ 2.718).
#' This is a symmetric p-value merger, admissible under general dependence by Vovk and Wang (2020).
#'
#' @param pvals A numeric vector of p-values in (0, 1].
#'
#' @return A numeric scalar: the combined p-value.
#' @references Vovk, V., & Wang, R. (2020). Combining p-values via averaging. \emph{Biometrika}, 107(4), 791–808.
#' @export
Geomean_pcombine <- function(pvals) {
  pvals <- pmin(pmax(pvals, .Machine$double.eps), 1)  # clamp to avoid log(0)
  gmean <- exp(mean(log(pvals)))
  pval <- min(exp(1) * gmean, 1)
  pval
}

#' Gmean P-Value Combination (Vovk et al., 2022)
#'
#' This function combines a vector of p-values using the Grid Gmean method 
#' described in Vovk, Wang, and Wang (2022, Annals of Statistics, Proposition 5).
#' It calculates multiple M-mean p-values over a geometric grid of powers and returns 
#' the minimum adjusted p-value.
#'
#' @param pvals Numeric vector of p-values in (0, 1].
#' @param r = 20 default
#'
#' @return A list with:
#' \describe{
#'   \item{min_p}{Combined p-value from the grid harmonic rule.}
#'   \item{r_grid}{Grid of r values used.}
#'   \item{F_grid}{F values for each r.}
#'   \item{selected_r}{Selected r.}
#' }
#' 
#' @references Vovk, V., Wang, R., & Wang, R. (2022). Admissible ways of merging
#' p-values under arbitrary dependence.
#'
#' @export
Genmean_pcombine <- function(pvals, r = 20) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n_p <- length(pvals)
  r = mean(1/(r-1),n_p-1)
  M_value = (mean(pvals^r))^(1/r)
  pval <- min(min(r+1,n_p)^(1/r) * M_value, 1)
  pval
}

#' Gmean P-Value Combination (Vovk et al., 2022)
#'
#' This function combines a vector of p-values using the Grid Gmean method 
#' described in Vovk, Wang, and Wang (2022, Annals of Statistics, Proposition 5).
#' It calculates multiple M-mean p-values over a geometric grid of powers and returns 
#' the minimum adjusted p-value.
#'
#' @param pvals Numeric vector of p-values in (0, 1].
#' @param r = 20 default
#'
#' @return A list with:
#' \describe{
#'   \item{min_p}{Combined p-value from the grid harmonic rule.}
#'   \item{r_grid}{Grid of r values used.}
#'   \item{F_grid}{F values for each r.}
#'   \item{selected_r}{Selected r.}
#' }
#' 
#' @references Vovk, V., Wang, R., & Wang, R. (2022). Admissible ways of merging
#' p-values under arbitrary dependence.
#'
#' @export
Genmean_rneg_pcombine <- function(pvals, r = -20) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  
  n_p <- length(pvals)
  M_value = (mean(pvals^r))^(1/r)
  pval <- min((r/(r+1))* n_p^(1+1/r) * M_value, 1)
  pval
}

#' IU P-value Combination Test (Spreng and Urga, 2022)
#'
#' Combines a vector of p-values.
#'
#' @param pvals A numeric vector of p-values (all in (0, 1]).
#' @param r Power parameter (default = 20).
#' @param alpha Significance level (default = 0.05).
#'
#' @return A list with test statistic, critical value, rejection indicator, and combined p-value.
#' @references Spreng, L., & Urga, G. (2022). Combining p-values for Multivariate Predictive Ability Testing.
#' @export
iu_pcombine <- function(pvals, r = 20) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n <- length(pvals)
  Prn <- (mean(pvals^(-r)))^(1 / r)
  pval <- min(r / ((r - 1) * Prn), 1)
  pval
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
bonferroni_pcombine <- function(pvals) {
  if (any(is.na(pvals)) || any(pvals < 0) || any(pvals > 1)) {
    stop("All p-values must be in [0, 1] and non-NA.")
  }
  n <- length(pvals)
  pval <- min(n * min(pvals), 1)
  pval
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
cauchy_pcombine <- function(p, weights = NULL) {
  if (any(p < 0 | p > 1)) stop("All p-values must lie strictly in [0,1].")
  
  n <- length(p)
  if (is.null(weights)) {
    weights <- rep(1 / n, n)
  } else {
    if (length(weights) != n || any(weights < 0)) stop("Weights must be non-negative and match p-values in length.")
    weights <- weights / sum(weights)  # normalize
  }

  T_stat <- sum(weights * tan((0.5 - p) * pi))
  pval <- min(max(0.5 - atan(T_stat) / pi, 0), 1)
  pval
}