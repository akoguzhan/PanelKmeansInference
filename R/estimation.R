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
    kappa = 1.5,
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

      gamma <- result$final_cluster
      theta <- result$centers[[result$iter]]
      gamma_tot <- gamma[match(id, unique(id))]
      D <- model.matrix(~ factor(gamma_tot) - 1)
      U <- Z - D %*% theta
      Sigma_hat <- t(U) %*% U / NT
      K_tot <- length(c(theta)) + length(gamma)

      # BIC <- log(det(Sigma_hat)) + (Ki * P + N) * kappa * sqrt(NT) / NT
      # BIC <- log(det(Sigma_hat)) + (Ki * P + N) * kappa * 1 / (5 * log(Tobs) * Tobs^(1/8)) # Liu et al. (2020)
      BIC <- log(det(Sigma_hat)) + (Ki * P + N) * kappa * log(NT) / (NT)
      # BIC <- log(det(Sigma_hat)) + (Ki * P + N) * kappa * log(N_eff * Tobs) / (N_eff * Tobs)
      # BIC <- log(det(Sigma_hat)) + (Ki * P + N) * 2 * log(log(NT)) / NT
      # BIC <- log(det(Sigma_hat)) + (Ki * P + N) * 2 / NT

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

estimate_alpha <- function(Z, id, time) {
  N <- length(unique(id))
  Tobs <- length(unique(time))

  Z_array <- array(NA, dim = c(N, Tobs, ncol(Z)))
  for (p in 1:ncol(Z)) {
    Z_array[, , p] <- tapply(Z[, p], list(id, time), mean)
  }

  alpha_vec <- sapply(1:ncol(Z), function(p) {
    cs_avg <- colMeans(Z_array[, , p])
    var_cs_avg <- var(cs_avg)
    alpha <- 1 + 0.5 * log(var_cs_avg) / log(N)
    return(alpha)
  })
  mean(alpha_vec, na.rm = TRUE)
}

panel_kmeans_cv_time <- function(
    Z, id, time,
    Kmax,
    nfold = 5,
    Ninit = 10, iter.max = 10,
    n_cores = 1) {
  
  if (!is.matrix(Z)) stop("Z must be a matrix")
  if (length(id) != nrow(Z)) stop("Length of id must match number of rows in Z")
  if (length(time) != nrow(Z)) stop("Length of time must match number of rows in Z")

  N <- length(unique(id))
  T_all <- sort(unique(time))
  Tobs <- length(T_all)
  P <- ncol(Z)

  # Split time into folds
  time_folds <- split(T_all, cut(seq_along(T_all), breaks = nfold, labels = FALSE))

  cv_errors <- numeric(Kmax - 1)
  names(cv_errors) <- paste0("K=", 2:Kmax)

  for (K in 2:Kmax) {
    fold_errors <- numeric(nfold)

    for (f in 1:nfold) {
      test_times <- time_folds[[f]]
      train_idx <- !(time %in% test_times)
      test_idx  <- (time %in% test_times)

      # Estimate K-means on training set
      res <- panel_kmeans_estimation(
        Z = Z[train_idx, , drop = FALSE],
        id = id[train_idx],
        time = time[train_idx],
        K = K,
        Ninit = Ninit,
        iter.max = iter.max,
        n_cores = n_cores
      )

      if (is.null(res)) {
        fold_errors[f] <- Inf
        next
      }

      # Get cluster centers and assignment
      gamma_train <- res$final_cluster
      theta_hat <- res$centers[[res$iter]]
      gamma_tot <- gamma_train[match(id[test_idx], sort(unique(id)))]

      # Compute test error
      D_test <- model.matrix(~ factor(gamma_tot) - 1)
      U_test <- Z[test_idx, , drop = FALSE] - D_test %*% theta_hat
      fold_errors[f] <- mean(rowSums(U_test^2))
    }

    cv_errors[K - 1] <- mean(fold_errors)
  }

  best_K <- which.min(cv_errors) + 1

  final_result <- panel_kmeans_estimation(
    Z = Z, id = id, time = time,
    K = best_K, Ninit = Ninit,
    iter.max = iter.max, n_cores = n_cores
  )

  final_result$CV_selected_K <- best_K
  final_result
}
