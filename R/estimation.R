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

  if (!is.null(K) && !is.null(Kmax)) {
    stop("Specify only one of K or Kmax, not both.")
  }
  if (!is.matrix(Z)) stop("Z must be a matrix")
  if (length(id) != nrow(Z)) stop("Length of id must match number of rows in Z")
  if (length(time) != nrow(Z)) stop("Length of time must match number of rows in Z")

  N <- length(unique(id))
  Tobs <- length(unique(time))
  P <- ncol(Z)
  NT <- N * Tobs
  kappa <- 1.5

  estimate_kmeans <- function(K) {

    # Parallel or sequential over Ninit initializations
    run_one <- function(init) {
      # Force every cluster to be represented at least once

      gamma <- rep(1:K, length.out = N)
      gamma <- sample(gamma, N)
      gamma_tot <- gamma[match(id, unique(id))]

      iter <- 1
      cluster_assign_list <- list()
      centroid_list <- list()
      objective_value <- list()
      converged <- FALSE

      while (iter <= iter.max && !converged) {
        # Compute cluster means
        theta <- aggregate_matrix(Z, gamma_tot)
        centroid_list[[iter]] <- theta

        # Assign clusters
        distance_matrix <- squared_distance(theta, Z)
        id_unique <- unique(id)
        objectives <- aggregate_matrix(t(distance_matrix), match(id, id_unique)) * Tobs
        gamma_new <- apply(objectives, 1, which.min)

        # Check if all clusters are present
        if (length(unique(gamma_new)) != K) {
          return(NULL) # Abandon this run
        }

        cluster_assign_list[[iter]] <- gamma_new
        objective_value[[iter]] <- sum(apply(objectives, 2, min))

        # Check for convergence
        if (iter > 1 && all(gamma_new == gamma)) {
          converged <- TRUE
        }
        gamma <- gamma_new
        gamma_tot <- gamma[match(id, unique(id))]
        iter <- iter + 1
      }

      gamma_final <- cluster_assign_list[[iter - 1]]
      names(gamma_final) <- as.character(seq_along(gamma_final))

      list(
        clusters = cluster_assign_list,
        centers = centroid_list,
        objectives = objective_value,
        iter = iter - 1,
        final_cluster = gamma_final
      )
    }

    if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
      if (.Platform$OS.type == "windows") {
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl))
        results <- parallel::parLapply(cl, 1:Ninit, run_one)
      } else {
        results <- parallel::mclapply(1:Ninit, run_one, mc.cores = n_cores)
      }
    } else {
      results <- lapply(1:Ninit, run_one)
    }
    # Remove failed runs
    results <- Filter(Negate(is.null), results)
    if (length(results) == 0) {
      warning(paste("All", Ninit, "initializations failed for K =", K))
      return(NULL)
    }
    objs <- sapply(results, function(res) min(unlist(res$objectives)))
    best_idx <- which.min(objs)
    return(results[[best_idx]])
  }
  
  if (!is.null(Kmax)) {
    BIC_opt <- Inf
    best_result <- NULL
    for (Ki in 2:Kmax) {
      result <- estimate_kmeans(Ki)
      if (is.null(result)) next

      theta <- result$centers[[result$iter]]
      gamma <- result$final_cluster
      gamma_tot <- gamma[match(id, unique(id))]
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
    result <- best_result
  } else if (!is.null(K)) {
    result <- estimate_kmeans(K)
  } else {
    stop("Either K or Kmax must be provided.")
  } 
  return(result)
}