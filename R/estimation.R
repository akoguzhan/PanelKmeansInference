#' Panel K-means Estimation with Optional BIC Selection (Parallelized)
#'
#' This function performs K-means clustering on NT x P matrix panel data. If `Kmax` is provided,
#' it performs model selection using BIC over 2 to `Kmax` clusters. Otherwise, clustering is
#' performed using the provided `K`. Parallelization is used for random initializations.
#'
#' @param Z NT x P matrix of variables for clustering (typically stacked by time)
#' @param id NT-length vector of unit identifiers (repeating across time)
#' @param time NT-length vector of time identifiers (repeating across units)
#' @param K Integer; number of clusters. Required unless Kmax is provided.
#' @param Kmax Optional; if provided, perform BIC-based selection from 2 to Kmax clusters.
#' @param Ninit Number of random initializations
#' @param iter.max Maximum number of Lloyd iterations
#' @param n_cores Number of cores for parallelization (default: 1, i.e., no parallel)
#'
#' @return A list with: clusters, centers, objectives, iter, final_cluster, BIC_selected_K (if Kmax is used)
#' @export
panel_kmeans_estimation <- function(
    Z, id, time,
    K = NULL, Kmax = NULL,
    Ninit = 10, iter.max = 10,
    n_cores = 1) {

  if (!is.matrix(Z)) stop("Z must be a matrix")
  if (length(id) != nrow(Z)) stop("Length of id must match number of rows in Z")
  if (length(time) != nrow(Z)) stop("Length of time must match number of rows in Z")

  N <- length(unique(id))
  Tobs <- length(unique(time))
  P <- ncol(Z)
  NT <- N * Tobs
  kappa <- 3

  estimate_kmeans <- function(K) {
    results <- list()
    reps <- 1
    while (reps <= Ninit) {
      cluster_assign_list <- vector("list", iter.max + 1)
      centroid_list <- vector("list", iter.max + 1)
      objective_value <- vector("list", iter.max + 1)

      iter <- 0
      gamma <- sample(1:K, size = N, replace = TRUE)
      gamma_tot <- rep(gamma, each = Tobs)
      theta <- aggregate_matrix(Z, gamma_tot)

      iter <- iter + 1
      distance_matrix <- squared_distance(theta, Z)
      id_unique <- unique(id)
      objectives <- aggregate_matrix(t(distance_matrix), match(id, id_unique)) * Tobs
      gamma <- apply(objectives, 1, which.min)
      gamma_tot <- rep(gamma, each = Tobs)

      if (length(unique(gamma)) != K) next

      centroid_list[[iter]] <- theta
      cluster_assign_list[[iter]] <- gamma
      objective_value[[iter]] <- sum(apply(objectives, 2, min))

      same_cluster <- FALSE
      while (iter <= iter.max && !same_cluster) {
        iter <- iter + 1
        theta <- aggregate_matrix(Z, gamma_tot)
        distance_matrix <- squared_distance(theta, Z)
        objectives <- aggregate_matrix(t(distance_matrix), match(id, id_unique)) * Tobs
        gamma <- apply(objectives, 1, which.min)
        gamma_tot <- rep(gamma, each = Tobs)

        centroid_list[[iter]] <- theta
        cluster_assign_list[[iter]] <- gamma
        same_cluster <- same_cl(gamma, cluster_assign_list[[iter - 1]], K)
        objective_value[[iter]] <- sum(apply(objectives, 2, min))
      }

      # Only keep results with exactly K clusters
      if (length(unique(gamma)) == K) {
        cluster_assign_list <- cluster_assign_list[1:iter]
        centroid_list <- centroid_list[1:iter]
        objective_value <- objective_value[1:iter]
        gamma_final <- cluster_assign_list[[iter]]
        names(gamma_final) <- as.character(seq_along(gamma_final))

        results[[length(results) + 1]] <- list(
          clusters = cluster_assign_list,
          centers = centroid_list,
          objectives = objective_value,
          iter = iter,
          final_cluster = gamma_final
        )
        reps <- reps + 1
      }
    }
    # Select the best result (lowest objective)
    objs <- sapply(results, function(res) min(unlist(res$objectives)))
    best_idx <- which.min(objs)
    return(results[[best_idx]])
  }

  if (!is.null(Kmax)) {
    BIC_opt <- Inf
    best_result <- NULL
    for (Ki in 2:Kmax) {
      result <- estimate_kmeans(Ki)
      theta <- result$centers[[result$iter]]
      gamma <- result$final_cluster
      gamma_tot <- rep(gamma, each = Tobs)
      D <- model.matrix(~ factor(gamma_tot) - 1)
      U <- Z - D %*% theta
      Sigma_hat <- t(U) %*% U / NT
      K_tot <- length(c(theta)) + length(gamma)
      BIC <- log(det(Sigma_hat)) + K_tot * kappa * log(NT) / NT
      if (BIC < BIC_opt) {
        BIC_opt <- BIC
        best_result <- result
        best_result$BIC_selected_K <- Ki
      }
    }
    return(best_result)
  }

  if (is.null(K)) stop("Either K or Kmax must be specified.")
  return(estimate_kmeans(K))
}