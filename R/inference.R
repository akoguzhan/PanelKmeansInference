#' Inference Test for Equality of Two Cluster Means in Panel K-means
#'
#' Tests whether two estimated cluster means are equal using standardized
#' differences and the truncated chi-squared distribution.
#' This version uses the \code{truncdist} package (with the gamma workaround)
#' to compute the p-value when the truncation interval is a single interval.
#'
#' @details
#' This function requires the \code{truncdist} package. The truncated chi-squared
#' distribution is handled via the gamma distribution with shape = df/2 and scale = 2.
#' If the truncation set is not a single interval, it falls back to the original method.
#'
#' @param Z NT x P matrix of data used in inference
#' @param id Numeric vector of unit identifiers (length NT)
#' @param time Numeric vector of time identifiers (length NT)
#' @param k1 First cluster index
#' @param k2 Second cluster index
#' @param required_quantities Optional list from required_quantities(); computed internally if NULL
#' @param estimated_k_means Output of panel_kmeans_estimation()
#'
#' @return List with final_interval (Intervals), test_stat (numeric), pval (numeric)
#' 
#' 
panel_kmeans_inference <- function(
    Z, id, time,
    required_quantities,
    estimated_k_means
) {

  Tobs <- required_quantities$Tobs
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
  pval <- TChisqRatioApprox(P, numer, denom)

  list(
    final_interval = S,
    test_stat = test_stat,
    pval = pval
  )
}

#' Selective Inference Test for Homogeneity of Cluster Means (Panel K-means)
#'
#' Performs a selective inference test for the global null hypothesis of cluster mean homogeneity
#' in panel K-means clustering. The function internally estimates clusters and the long-run variance,
#' then combines pairwise selective p-values using a specified combination function.
#'
#' @param Z NT x P matrix of data used in inference
#' @param id Numeric vector of unit identifiers (length NT)
#' @param time Numeric vector of time identifiers (length NT)
#' @param K Integer; number of clusters. Required unless Kmax is provided.
#' @param Kmax Optional; if provided, perform BIC-based selection from 2 to Kmax clusters.
#' @param Ninit Number of random initializations for k-means (default: 10)
#' @param iter.max Maximum number of Lloyd iterations (default: 10)
#' @param lrv String; long-run variance estimator to use (e.g., "EWC")
#' @param lrv_par Optional; parameters for the long-run variance estimator
#' @param pairs Optional matrix of cluster index pairs to test (each row: c(k1, k2)). If NULL, all unique pairs are tested.
#' @param pcombine_fun Function to combine p-values: "pmean", "porder", "pSimes", "pharmonic", "pCauchy"
#' @param r Used by pmean
#'
#' @return List with pairwise_pvalues, pairs, pvalue_combination
#' @export
panel_homogeneity_test <- function(
    Z, id, time,
    K = NULL,
    Kmax = NULL,
    Ninit = 10,
    iter.max = 10,
    lrv = "EWC",
    lrv_par = NULL,
    pairs = NULL,
    pcombine_fun = "Genmean_neq",
    r = -20,
    n_cores = 1
) {
  # Estimate clusters
  estimated_k_means <- panel_kmeans_estimation(
    Z = Z,
    id = id,
    time = time,
    K = K,
    Kmax = Kmax,
    Ninit = Ninit,
    iter.max = iter.max,
    n_cores = n_cores
  )
  K_est <- length(unique(estimated_k_means$final_cluster))
  gamma <- estimated_k_means$final_cluster

  K_used <- if (!is.null(K)) K else estimated_k_means$BIC_selected_K
  N <- length(unique(id))
  Tobs <- length(unique(time))
  P <- ncol(Z)
  KP <- K_used * P

  # 2. Calculate long-run variance matrix (OmegaHat) for the clustered means
  id_unique <- sort(unique(id))
  gamma_long <- gamma[match(id, id_unique)]
  Z_by_group <- matrix(NA, nrow = Tobs, ncol = KP)
  for (k in 1:K_used) {
    idx_k <- gamma_long == k
    Z_k <- Z[idx_k, , drop = FALSE]
    time_k <- time[idx_k]
    Zbar_k <- aggregate_matrix(Z_k, time_k)  # Tobs x P
    Z_by_group[, ((k - 1) * P + 1):(k * P)] <- Zbar_k
  }
  if (lrv == "NeweyWest") {
    OmegaHat <- NeweyWest(Z_by_group, lrv_par = lrv_par)
  } else if (lrv == "EWC") {
    OmegaHat <- EWC(Z_by_group, lrv_par = lrv_par)$S
  } else {
    stop("Unknown lrv. Use 'EWC' or 'NeweyWest'.")
  }

  # Determine which pairs to test
  if (is.null(pairs)) {
    pair_idx <- t(combn(K_est, 2))
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
      required_quantities = required_quantities(Z, id, time, k1, k2, gamma, OmegaHat),
      estimated_k_means = estimated_k_means
    )$pval
    pairs_out[[j]] <- c(k1, k2)
  }

  # Use your cauchy_combine function if requested
  merge_result <- switch(
    pcombine_fun,
    Geomean     = Geomean_pcombine(pvalues),
    Genmean     = Genmean_pcombine(pvalues, r = r),
    Genmean_neq = Genmean_rneg_pcombine(pvalues, r = r),
    iu          = iu_pcombine(pvalues, r = r),
    bonferroni  = bonferroni_pcombine(pvalues),
    cauchy      = cauchy_pcombine(pvalues),
    stop("Invalid pcombine_fun specified.")
  )

  list(
    pairwise_pvalues = pvalues,
    pairs = pairs_out,
    pval = merge_result,
    clustering = gamma
  )
}

#' Selective Inference: Test if Each Cluster Mean is Zero
#'
#' Performs a selective inference test for H0: mu_k = 0 for each cluster k.
#'
#' @param Z NT x P matrix of data.
#' @param id Vector of unit identifiers (length NT).
#' @param time Vector of time identifiers (length NT).
#' @param K Integer; number of clusters. Required unless Kmax is provided.
#' @param Kmax Optional; if provided, perform BIC-based selection from 2 to Kmax clusters.
#' @param Ninit Number of k-means initializations.
#' @param iter.max Max iterations for k-means.
#' @param lrv Long-run variance method: \"EWC\" or \"NeweyWest\".
#' @param lrv_par Parameter for long-run variance estimator.
#' @param n_cores Number of cores to use.
#'
#' @return List with individual p-values and estimated clustering.
#' @export
panel_bycluster_tests <- function(
  Z, id, time,
  K = NULL,
  Kmax = NULL,
  Ninit = 10,
  iter.max = 10,
  lrv = "EWC",
  lrv_par = NULL,
  n_cores = 1
) {
  # Step 1: Estimate clusters
  bycluster_res <- panel_kmeans_estimation(
    Z = Z,
    id = id,
    time = time,
    K = K,
    Kmax = Kmax,
    Ninit = Ninit,
    iter.max = iter.max,
    n_cores = n_cores
  )
  gamma <- bycluster_res$final_cluster
  K_est <- length(unique(gamma))

  N <- length(unique(id))
  Tobs <- length(unique(time))
  P <- ncol(Z)
  KP <- K_est * P

  # Step 2: Long-run variance matrix for group means
  id_unique <- sort(unique(id))
  gamma_long <- gamma[match(id, id_unique)]
  Z_by_group <- matrix(NA, nrow = Tobs, ncol = KP)
  for (k in 1:K_est) {
    idx_k <- gamma_long == k
    Z_k <- Z[idx_k, , drop = FALSE]
    time_k <- time[idx_k]
    Zbar_k <- aggregate_matrix(Z_k, time_k)  # Tobs x P
    Z_by_group[, ((k - 1) * P + 1):(k * P)] <- Zbar_k
  }

  OmegaHat <- switch(
    lrv,
    EWC = EWC(Z_by_group, lrv_par = lrv_par)$S,
    NeweyWest = NeweyWest(Z_by_group, lrv_par = lrv_par),
    stop("Unknown lrv. Use 'EWC' or 'NeweyWest'.")
  )

  # Step 3: Loop over clusters, test H0: mu_k = 0
  pvals <- numeric(K_est)
  for (k in 1:K_est) {
    rq <- required_quantities(
      Z = Z, id = id, time = time,
      k1 = k, k2 = 0,
      gamma = gamma,
      OmegaHat = OmegaHat
    )

    pvals[k] <- panel_kmeans_inference(
      Z = Z, id = id, time = time,
      required_quantities = rq,
      estimated_k_means = bycluster_res
    )$pval
  }

  return(list(
    bycluster_pvals = pvals,
    clustering = gamma
  ))
}