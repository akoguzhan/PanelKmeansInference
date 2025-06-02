#' Inference Test for Equality of Two Cluster Means in Panel K-means
#'
#' Tests whether two estimated cluster means are equal using standardized
#' differences and truncated chi-squared distribution.
#'
#' @param Z NT x P matrix of data used in inference
#' @param id Numeric vector of unit identifiers (length NT)
#' @param time Numeric vector of time identifiers (length NT)
#' @param k1 First cluster index
#' @param k2 Second cluster index
#' @param OmegaHat GK x GK long-run variance matrix
#' @param required_quantities Optional list from required_quantities(); computed internally if NULL
#' @param estimated_k_means Output of panel_kmeans_estimation()
#'
#' @return List with final_interval (Intervals), test_stat (numeric), pval (numeric)
#'
panel_kmeans_inference <- function(
    Z, id, time,
    k1, k2,
    OmegaHat,
    required_quantities = NULL,
    estimated_k_means
) {
  if (is.null(required_quantities)) {
    required_quantities <- required_quantities(
      Z = Z,
      id = id,
      time = time,
      k1 = k1,
      k2 = k2,
      gamma = estimated_k_means$final_cluster,
      OmegaHat = OmegaHat
    )
  }

  Tobs <- required_quantities$Tobs
  K <- required_quantities$K
  P <- required_quantities$P
  Omega_ZTv_norm <- required_quantities$Omega_ZTv_norm

  test_stat <- sqrt(Tobs) * Omega_ZTv_norm
  gestat <- intervals::Intervals(c(test_stat^2, Inf))

  S <- compute_S(required_quantities, estimated_k_means)
  S <- suppressWarnings(intervals::interval_intersection(S, intervals::Intervals(c(0, Inf))))

  if (nrow(S) > 1) {
    checks <- apply(S, 1, function(row) {
      nrow(suppressWarnings(intervals::interval_intersection(intervals::Intervals(row), gestat))) > 0
    })
    if (any(checks)) {
      S <- S[checks, , drop = FALSE]
    }
  }

  denom <- S^2
  numer <- suppressWarnings(intervals::interval_intersection(gestat, denom))
  pval <- TChisqRatioApprox(K * P, numer, denom)

  list(
    final_interval = S,
    test_stat = test_stat,
    pval = pval
  )
}

#' Test Homogeneity of Cluster Means Using Combined Pairwise Inference
#'
#' Combines pairwise tests of cluster mean equality using pmerge package.
#'
#' @param Z NT x P matrix
#' @param id Numeric vector of unit identifiers
#' @param time Numeric vector of time identifiers
#' @param OmegaHat Long-run variance matrix (GK x GK)
#' @param estimated_k_means Output from panel_kmeans_estimation()
#' @param pairs Optional matrix of cluster index pairs to test (each row: c(k1, k2)). If NULL, all unique pairs are tested.
#' @param pcombine_fun Function to combine p-values: "pmean", "porder", "pSimes", "pharmonic", "pCauchy"
#' @param method Variant used by pmerge (e.g., "H1", "I", "G", etc.)
#' @param order_k Used by porder
#' @param r Used by pmean
#' @param epi Used by pCauchy
#'
#' @return List with combined_stat, pairwise_pvalues, pairs, pvalue_combination
#' @export
#'
panel_homogeneity_test <- function(
    Z, id, time,
    OmegaHat,
    estimated_k_means,
    pairs = NULL,
    pcombine_fun = "pharmonic",
    method = "H2",
    order_k = NULL,
    r = NULL,
    epi = NULL
) {
  K <- length(unique(estimated_k_means$final_cluster))
  # Determine which pairs to test
  if (is.null(pairs)) {
    pair_idx <- t(combn(K, 2))
  } else {
    pair_idx <- as.matrix(pairs)
    if (ncol(pair_idx) != 2) stop("pairs must be a matrix with two columns (k1, k2).")
  }
  npairs <- nrow(pair_idx)
  pvalues <- numeric(npairs)
  pairs_out <- vector("list", npairs)

  for (j in seq_len(npairs)) {
    k1 <- pair_idx[j, 1]
    k2 <- pair_idx[j, 2]
    pvalues[j] <- panel_kmeans_inference(
      Z = Z,
      id = id,
      time = time,
      k1 = k1,
      k2 = k2,
      OmegaHat = OmegaHat,
      estimated_k_means = estimated_k_means
    )$pval
    pairs_out[[j]] <- c(k1, k2)
  }

  merge_result <- switch(
    pcombine_fun,
    pmean     = pmerge::pmean(p = pvalues, r = r, dependence = method),
    porder    = pmerge::porder(p = pvalues, k = order_k),
    pSimes    = pmerge::pSimes(p = pvalues, method = method),
    pharmonic = pmerge::pharmonic(p = pvalues, method = method),
    stop("Invalid pcombine_fun specified.")
  )

  list(
    pairwise_pvalues = pvalues,
    pairs = pairs_out,
    pvalue_combination = merge_result
  )
}
