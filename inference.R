#' Selective Inference for Panel K-means Cluster Means
#'
#' This function conducts a selective inference test for cluster-specific means 
#' estimated via panel K-means, adjusting for the cluster selection process.
#' If \code{g2} is omitted, the null hypothesis is that the mean of cluster \code{g1}
#' is zero.
#'
#' @param df A balanced panel data.frame
#' @param id Character; name of the unit identifier column
#' @param time Character; name of the time identifier column
#' @param Z_names Character vector of variables used for clustering/inference
#' @param K Integer; number of clusters
#' @param g1 Integer; index of the first cluster
#' @param g2 Integer (optional); index of second cluster. If \code{NULL}, test is against zero.
#' @param lr_method Character; either "EWC" (default) or "NeweyWest" for long-run variance
#' @param lr_param Integer; smoothing parameter (B for EWC, maxlag for NeweyWest). Optional
#' @param Ninit Integer; number of K-means initializations
#' @param iter.max Integer; max iterations in K-means
#'
#' @return A list with:
#' \describe{
#'   \item{\code{test_stat}}{Standardized test statistic}
#'   \item{\code{final_interval}}{Truncation set (Intervals object)}
#'   \item{\code{pval}}{Selective p-value}
#'   \item{\code{clusters}}{Final cluster assignments}
#' }
#' @export
cluster_test <- function(df, id, time, Z_names,
                         K, g1, g2 = NULL,
                         lr_method = "EWC", lr_param = NULL,
                         Ninit = 10, iter.max = 10) {
  # Step 1: Convert to matrices
  panel_list <- panel_data_to_matrix(df, id, time, Z_names)
  Z <- panel_list$Z
  id_vec <- panel_list$id
  time_vec <- panel_list$time
  
  # Step 2: Estimate panel K-means
  km_out <- panel_kmeans_estimation(Z = Z, id = id_vec,
                                    K = K, Ninit = Ninit, iter.max = iter.max)
  
  # Step 3: Compute long-run variance (GK x GK)
  Tobs <- length(unique(time_vec))
  P <- length(Z_names)
  Z_df <- as.data.frame(Z)
  Z_df$time <- time_vec
  Z_df$cluster <- rep(km_out$final_cluster, each = Tobs)
  
  Zbar_df <- aggregate(Z_df[, Z_names], 
                       by = list(time = Z_df$time, cluster = Z_df$cluster), 
                       FUN = mean)
  Zbar_mat <- as.matrix(reshape(Zbar_df, 
                                idvar = "time", 
                                timevar = "cluster", 
                                direction = "wide")[, -1])
  
  if (lr_method == "EWC") {
    if (is.null(lr_param)) lr_param <- Tobs
    OmegaHat <- EWC(Zbar_mat, B = lr_param)
  } else if (lr_method == "NeweyWest") {
    if (is.null(lr_param)) lr_param <- floor(1.3 * sqrt(Tobs))
    OmegaHat <- NeweyWest(Zbar_mat, maxlag = lr_param)
  } else {
    stop("Invalid lr_method. Use 'EWC' or 'NeweyWest'.")
  }
  
  # Step 4: Compute required quantities
  rq <- required_quantities(Z = Z,
                            id = id_vec, time = time_vec,
                            g1 = g1, g2 = g2,
                            gamma = km_out$final_cluster,
                            OmegaHat = OmegaHat)
  
  # Step 5: Inference test
  out <- panel_kmeans_inference(Z = Z,
                                id = id_vec, time = time_vec,
                                g1 = g1, g2 = g2,
                                OmegaHat = OmegaHat,
                                estimated_k_means = km_out,
                                required_quantities = rq)
  
  return(list(
    test_stat = out$test_stat,
    final_interval = out$final_interval,
    pval = out$pval,
    clusters = km_out$final_cluster
  ))
}

#' Inference Test for Equality of Two Cluster Means in Panel K-means
#'
#' Tests whether two estimated cluster means are equal using standardized
#' differences and truncated chi-squared distribution.
#'
#' @param Z NT x P matrix of data used in inference
#' @param id Numeric vector of unit identifiers (length NT)
#' @param g1 First cluster index
#' @param g2 Second cluster index
#' @param OmegaHat GK x GK long-run variance matrix
#' @param required_quantities Optional list from required_quantities(); computed internally if NULL
#' @param estimated_k_means Output of panel_kmeans_estimation()
#'
#' @return List with final_interval (Intervals), test_stat (numeric), pval (numeric)
#' @export
panel_kmeans_inference <- function(
    Z, id,
    g1, g2,
    OmegaHat,
    required_quantities = NULL,
    estimated_k_means
) {
  if (is.null(required_quantities)) {
    required_quantities <- required_quantities(
      Z = Z,
      id = id,
      g1 = g1,
      g2 = g2,
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
#' @param OmegaHat Long-run variance matrix (GK x GK)
#' @param estimated_k_means Output from panel_kmeans_estimation()
#' @param pcombine_fun Function to combine p-values: "pmean", "porder", "pSimes", "pharmonic", "pCauchy"
#' @param method Variant used by pmerge (e.g., "H1", "I", "G", etc.)
#' @param k Used by porder
#' @param r Used by pmean
#' @param epi Used by pCauchy
#'
#' @return List with combined_stat, pairwise_pvalues, pvalue_combination
#' @export
panel_homogeneity_test <- function(
    Z, id,
    OmegaHat,
    estimated_k_means,
    pcombine_fun = "pharmonic",
    method = "H1",
    k = NULL,
    r = NULL,
    epi = NULL
) {
  K <- length(unique(estimated_k_means$final_cluster))
  pvalues <- numeric(K - 1)
  
  for (g in 2:K) {
    pvalues[g - 1] <- panel_kmeans_inference(
      Z = Z,
      id = id,
      g1 = 1,
      g2 = g,
      OmegaHat = OmegaHat,
      estimated_k_means = estimated_k_means
    )$pval
  }
  
  merge_result <- switch(
    pcombine_fun,
    pmean     = pmerge::pmean(p = pvalues, r = r, dependence = method),
    porder    = pmerge::porder(p = pvalues, k = k),
    pSimes    = pmerge::pSimes(p = pvalues, method = method),
    pharmonic = pmerge::pharmonic(p = pvalues, method = method),
    pCauchy   = pmerge::pCauchy(p = pvalues, method = method, epi = epi),
    stop("Invalid pcombine_fun specified.")
  )
  
  list(
    combined_stat = merge_result$stat,
    pairwise_pvalues = pvalues,
    pvalue_combination = merge_result$pvalue
  )
}