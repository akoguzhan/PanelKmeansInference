#' Prepare Quantities for Truncation Set in Two-Group Inference
#'
#' This internal function computes key quantities used for selective inference in 
#' comparing cluster means. Inputs are assumed to be already in matrix form with
#' numeric `id` and `time` vectors of length NT.
#'
#' @param Z NT x P matrix of stacked panel data (P can be 1)
#' @param id Vector of unit identifiers (length NT, integers)
#' @param time Vector of time identifiers (length NT, integers)
#' @param k1 First cluster index (scalar)
#' @param k2 Second cluster index (scalar)
#' @param gamma Vector of length N with estimated cluster labels
#' @param OmegaHat GK x GK long-run variance-covariance matrix
#'
#' @return A list of required quantities for computing truncation intervals
#' @keywords internal
#'
required_quantities <- function(Z, id, time, k1, k2, gamma, OmegaHat) {
  # Ensure Z is a matrix (even if P = 1)
  Z <- as.matrix(Z)
  P <- ncol(Z)
  N <- length(unique(id))
  Tobs <- length(unique(time))
  K <- length(unique(gamma))
  NT <- N * Tobs

  # Expand cluster labels to NT
  gamma_tot <- rep(gamma, each = Tobs)

  # Weight vector v
  v <- rep(0, NT)
  v[gamma_tot == k1] <- 1 / sum(gamma_tot == k1)
  if (k2 != 0) {
    v[gamma_tot == k2] <- -1 / sum(gamma_tot == k2)
  }

  # Variance of the difference
  if (k2 == 0) {
    Vdiff <- OmegaHat[((k1 - 1) * P + 1):(k1 * P), ((k1 - 1) * P + 1):(k1 * P), drop = FALSE]
  } else {
    Vg1 <- OmegaHat[((k1 - 1) * P + 1):(k1 * P), ((k1 - 1) * P + 1):(k1 * P), drop = FALSE]
    Vg2 <- OmegaHat[((k2 - 1) * P + 1):(k2 * P), ((k2 - 1) * P + 1):(k2 * P), drop = FALSE]
    Cov12 <- OmegaHat[((k1 - 1) * P + 1):(k1 * P), ((k2 - 1) * P + 1):(k2 * P), drop = FALSE]
    Vdiff <- Vg1 + Vg2 - 2 * Cov12
  }

  # Statistics
  v_norm <- norm_vec(v)
  ZTv <- t(Z) %*% v
  ZTv_norm <- norm_vec(ZTv)
  dir_ZTv <- dir_vec(ZTv)
  Omega_ZTv_norm <- sqrt(as.numeric(t(ZTv) %*% solve(Vdiff) %*% ZTv))
  nr <- as.vector(ZTv_norm / Omega_ZTv_norm)

  return(list(
    id = id,
    time = time,
    P = P,
    K = K,
    N = N,
    Tobs = Tobs,
    Z = Z,
    v = v,
    v_norm = v_norm,
    ZTv = ZTv,
    ZTv_norm = ZTv_norm,
    dir_ZTv = dir_ZTv,
    Omega_ZTv_norm = Omega_ZTv_norm,
    nr = nr
  ))
}

#' Represent ||[z(ϕ)]_{it} - [z(ϕ)]_{jt}||_2^2 as a Quadratic in ϕ
#' In the current implementation, this function is not used (random initialization).
#'
#' Computes the coefficients of a quadratic function in φ representing the 
#' squared difference in transformed signals for units i and j.
#'
#' @param rq List of required quantities from \code{required_quantities()} (matrix-based)
#' @param i Integer; first unit index
#' @param j Integer; second unit index
#'
#' @return List with components: \code{quad}, \code{linear}, \code{constant}
#' @keywords internal
#'
norm_sq_phi_panel <- function(rq, i, j) {
  requireNamespace(Rfast)
  
  # Extract inputs
  id <- rq$id
  time <- rq$time
  Tobs <- rq$Tobs
  Z <- rq$Z
  v <- rq$v
  v_norm <- rq$v_norm
  ZTv <- rq$ZTv
  ZTv_norm <- rq$ZTv_norm
  dir_ZTv <- rq$dir_ZTv
  nr <- rq$nr
  
  # Extract unit-level data
  Z_i <- Z[id == i, , drop = FALSE]           # T x P matrix
  Z_j_bar <- colMeans(Z[id == j, , drop = FALSE])
  deltahat_i <- v[id == i & time == time[1]]
  deltahat_j <- v[id == j & time == time[1]]
  
  # Differences
  D_deltahat <- deltahat_i - deltahat_j
  D_Z <- Rfast::eachrow(Z_i, Z_j_bar, "-")    # T x P matrix
  dot_product <- rowSums(Rfast::eachrow(D_Z, t(dir_ZTv), "*"))
  
  # Quadratic coefficients
  a_ij <- nr^2 * (D_deltahat / (sqrt(Tobs) * v_norm^2))^2
  b_ijt <- 2 * nr * (
    (D_deltahat / (sqrt(Tobs) * v_norm^2)) * dot_product -
      (D_deltahat^2 / (sqrt(Tobs) * v_norm^4)) * ZTv_norm
  )
  D_resid <- Rfast::eachrow(D_Z, D_deltahat * ZTv / v_norm^2, "-")
  c_ijt <- rowSums(D_resid * D_resid)
  
  # Return quadratic representation
  list(
    quad = Tobs * a_ij,
    linear = sum(b_ijt),
    constant = sum(c_ijt)
  )
}

#' Quadratic Form of ||z_it(ϕ) - ClusterMean_k(ϕ)||² for Canonical Form
#'
#' Represent
#' \deqn{||[z(ϕ)]_{it} - T^{-1} ∑_{t} ∑_{j} w_j^{(m-1)}(k) [z(ϕ)]_{jt}||²}
#' as a quadratic function in ϕ.
#'
#' @param rq List from \code{required_quantities()} (matrix-based)
#' @param i Integer; unit index
#' @param k Integer; cluster index
#' @param last_centroids K x P matrix of centroids from previous k-means step
#' @param weighted_deltahat K-length vector of weighted δₖ
#'
#' @return A list with \code{quad}, \code{linear}, \code{constant} terms
#' @keywords internal
#'
norm_phi_canonical_panel <- function(rq, i, k,
                                     last_centroids,
                                     weighted_deltahat) {
  requireNamespace(Rfast)
  
  # Extract inputs
  id <- rq$id
  time <- rq$time
  Tobs <- rq$Tobs
  Z <- rq$Z
  v <- rq$v
  v_norm <- rq$v_norm
  ZTv <- rq$ZTv
  ZTv_norm <- rq$ZTv_norm
  dir_ZTv <- rq$dir_ZTv
  nr <- rq$nr
  
  # Observations for unit i
  Z_i <- Z[id == i, , drop = FALSE]
  deltahat_i <- v[id == i & time == time[1]]
  
  # Differences
  D_deltahat <- deltahat_i - weighted_deltahat[k]
  D_Z <- Rfast::eachrow(Z_i, last_centroids[k, ], "-")
  dot_product <- rowSums(Rfast::eachrow(D_Z, t(dir_ZTv), "*"))
  
  # Coefficients
  a_tilde <- nr^2 * (D_deltahat / (sqrt(Tobs) * v_norm^2))^2
  b_tilde <- 2 * nr * (
    (D_deltahat / (sqrt(Tobs) * v_norm^2)) * dot_product -
      (D_deltahat^2 / (sqrt(Tobs) * v_norm^4)) * ZTv_norm
  )
  resid <- Rfast::eachrow(D_Z, D_deltahat * ZTv / v_norm^2, "-")
  c_tilde <- rowSums(resid * resid)
  
  return(list(
    quad = Tobs * a_tilde,
    linear = sum(b_tilde),
    constant = sum(c_tilde)
  ))
}

#' Compute Cluster-Weighted δ̂ for Canonical Projection Step
#'
#' Computes \eqn{\sum_{j=1}^N w_j^{(m-1)}(k) \hat{\delta}_j} for each cluster k.
#'
#' @param rq Output of \code{required_quantities()} (matrix-based input)
#' @param last_cl Integer vector of length N; previous cluster assignments
#'
#' @return Numeric vector of length K with weighted δ̂ for each cluster
#' @keywords internal
#'
weighted_deltahat <- function(rq, last_cl) {
  K <- rq$K
  v <- rq$v
  time <- rq$time

  # Use time==first_period to subset one copy of δ̂ per unit
  v1 <- v[time == time[1]]

  weighted <- numeric(K)
  for (k in seq_len(K)) {
    indices <- which(last_cl == k)
    weighted[k] <- if (length(indices) > 0) {
      mean(v1[indices])
    } else {
      0
    }
  }

  return(weighted)
}

#' Compute Truncation Set S for Selective Inference
#'
#' @param rq List of required quantities computed from \code{required_quantities()}
#' @param estimated_k_means Result of \code{panel_kmeans_estimation()}
#'
#' @return An "Intervals" object representing the truncation set S
#' @keywords internal
#' 
compute_S <- function(rq, estimated_k_means) {
  requireNamespace("intervals")

  # Extract cluster assignment and center history
  cluster_path <- do.call(rbind, estimated_k_means$clusters)
  centroid_path <- estimated_k_means$centers
  M <- nrow(cluster_path) # number of assignment steps (including initialization)

  # Panel dimensions and identifiers
  N <- rq$N
  K <- rq$K
  ids <- if (!is.null(rq$ids)) rq$ids else unique(rq$id)
  Tobs <- rq$Tobs

  all_interval_lists <- list()

  # --- Initialization step: assignments given initial centers ---
  # The first row of cluster_path is the initial assignment (random), 
  # the first element of centroid_path is the initial centers.
  # The *second* row of cluster_path is the first Lloyd assignment.
  if (M >= 2) {
    initial_centers <- centroid_path[[1]]
    first_assignment <- cluster_path[2, ] # after first Lloyd update

    for (i in seq_len(N)) {
      id_i <- ids[i]
      assigned_k <- first_assignment[i]
      # For each possible cluster, encode the assignment event
      for (k in seq_len(K)) {
        # Only need to encode for k != assigned_k
        if (k != assigned_k) {
          # Compute quadratic for assigned cluster and for alternative cluster
          g_star_quad <- norm_phi_canonical_panel(rq, id_i, assigned_k, initial_centers, weighted_deltahat = rep(0, K))
          k_quad <- norm_phi_canonical_panel(rq, id_i, k, initial_centers, weighted_deltahat = rep(0, K))
          coeffs <- minus_quad_ineq(g_star_quad, k_quad)
          interval <- solve_one_ineq_complement(coeffs$quad, coeffs$linear, coeffs$constant)
          all_interval_lists[[length(all_interval_lists) + 1]] <- interval
        }
      }
    }
  }

  # --- Lloyd iterations: assignments given previous centers ---
  # For m = 3,...,M (since m=2 is the first Lloyd assignment)
  if (M > 2) {
    for (m in 3:M) {
      current_assignment <- cluster_path[m, ]
      previous_assignment <- cluster_path[m - 1, ]
      last_centroids <- centroid_path[[m - 1]]
      wdelta <- weighted_deltahat(rq, previous_assignment)

      for (i in seq_len(N)) {
        id_i <- ids[i]
        assigned_k <- current_assignment[i]
        for (k in seq_len(K)) {
          if (k != assigned_k) {
            g_star_quad <- norm_phi_canonical_panel(rq, id_i, assigned_k, last_centroids, wdelta)
            k_quad <- norm_phi_canonical_panel(rq, id_i, k, last_centroids, wdelta)
            coeffs <- minus_quad_ineq(g_star_quad, k_quad)
            interval <- solve_one_ineq_complement(coeffs$quad, coeffs$linear, coeffs$constant)
            all_interval_lists[[length(all_interval_lists) + 1]] <- interval
          }
        }
      }
    }
  }

  # --- Final intersection: compute S ---
  if (length(all_interval_lists) == 0) {
    stop("No truncation intervals were generated.")
  }
  raw_matrix <- matrix(do.call(c, all_interval_lists), ncol = 2, byrow = TRUE)
  complement_intervals <- intervals::reduce(intervals::Intervals(raw_matrix), check_valid = FALSE)
  final_S <- intervals::interval_complement(complement_intervals)

  return(final_S)
}

#' Implement the minus operation for two quadratic inequalities
#'
#' @keywords internal
#' 
#' @param quad1, first inequality
#' @param quad2, second inequality
#' 
minus_quad_ineq = function(quad1,quad2) {
  coef_list = list("quad" = quad1$quad-quad2$quad,
                   "linear" = quad1$linear-quad2$linear,
                   "constant" = quad1$constant-quad2$constant)
  return(coef_list)
}

#' Solve the roots of quadratic polynomials related to testing for a difference 
#' in means
#'
#' Solves \eqn{ax^2 + bx + c \ge 0}, then returns the complement of the solution 
#' set wrt to the real line, unless the complement is empty, in which case the 
#' function returns NA.
#'
#' @keywords internal
#'
#' @param A, B, C the coefficients of the quadratic equation.
#' @param tol if \eqn{|a|}, \eqn{|b|}, or \eqn{|c|} is not larger than tol, then 
#' treat it as zero.
#'
#' @return Returns an "Intervals" object containing NA or the complement of the 
#' solution set.
#' 
solve_one_ineq_complement = function(A,B,C,tol=1e-10) {
  # Computes the complement of the set {phi: B*phi + C <=  0},
  compute_linear_ineq_complement = function(B,C,tol=1e-8) {
    #  If B = 0
    if(abs(B) <= tol) {
      if(C <= tol) { # C <= 0: inequality is always satisfied
        return(c(0,0)) # all of real line
      } else { # C > 0: something has gone wrong -- no solution works
        warning("B = 0 and C > 0: B*phi + C <= 0 has no solution")
        return(c(-Inf,Inf)) # do not return any value
      }
    }
    
    # If B \neq 0
    ratio = -C/B
    # If B > 0:
    if(B > tol) {
      return(c(ratio,Inf))
    }
    if(B < tol) {
      return(c(-Inf, ratio))
    }
  }
  
  # A = 0?
  if(abs(A) <= tol) {
    return(compute_linear_ineq_complement(B,C,tol))
  }
  
  # We know A \neq 0
  discrim = B^2 - 4*A*C
  
  # If discriminant is small, we assume there is no root
  if(discrim <= tol) {
    if(A > tol) { # Parabola opens up: there is no solution
      return(c(-Inf,Inf))
    } else { # Parabola opens down: every x is a solution
      return(c(0, 0))
    }
  }
  
  # We now know that A =/= 0, and that there are two roots
  # we compute the roots using the suggestion outlined at
  # https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
  # for numerical stability purposes
  sqrt_discrim = sqrt(discrim)
  if (B >= tol){
    root_1 = (-B-sqrt_discrim)/(2*A)
    root_2 = (2*C)/(-B-sqrt_discrim)
    roots = sort(c(root_1, root_2))
  }else{
    root_1 = (-B+sqrt_discrim)/(2*A)
    root_2 = (2*C)/(-B+sqrt_discrim)
    roots = sort(c(root_1, root_2))
  }
  
  if(A > tol) {
    if(roots[1] > tol) {
      return(c(0, roots[1], roots[2], Inf))
    }
    
    if(roots[2] <= tol) {
      warning("something wrong with the discriminant calculation!")
      return(c(0,Inf))
    }
    
    return(c(roots[2], Inf))
  }
  
  # We now know that there are two roots, and parabola opens down (A < 0)
  if(roots[2] < -tol) {
    return(c(0, 0))
  }
  
  if(roots[1] > tol) {
    return(c(roots[1], roots[2]))
  }
  
  return(c(-Inf, roots[2]))
}

