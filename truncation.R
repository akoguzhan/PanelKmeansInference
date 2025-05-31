#' Prepare Quantities for Truncation Set in Two-Group Inference
#'
#' This internal function computes key quantities used for selective inference in 
#' comparing cluster means. Inputs are assumed to be already in matrix form with
#' numeric `id` and `time` vectors of length NT.
#'
#' @param Z NT x P matrix of stacked panel data
#' @param id Vector of unit identifiers (length NT, integers)
#' @param time Vector of time identifiers (length NT, integers)
#' @param g1 First cluster index (scalar)
#' @param g2 Second cluster index (scalar)
#' @param gamma Vector of length N with estimated cluster labels
#' @param OmegaHat GK x GK long-run variance-covariance matrix
#'
#' @return A list of required quantities for computing truncation intervals
#' @keywords internal
#'
required_quantities <- function(Z, id, time, g1, g2, gamma, OmegaHat) {
  # Scalars
  P <- ncol(Z)
  N <- length(unique(id))
  Tobs <- length(unique(time))
  K <- length(unique(gamma))
  NT <- N * Tobs
  
  # Expand cluster labels to NT
  gamma_tot <- rep(gamma, each = Tobs)
  
  # Weight vector v
  v <- rep(0, NT)
  v[gamma_tot == g1] <- 1 / sum(gamma_tot == g1)
  if (g2 != 0) {
    v[gamma_tot == g2] <- -1 / sum(gamma_tot == g2)
  }
  
  # Variance of the difference
  if (g2 == 0) {
    Vdiff <- OmegaHat[((g1 - 1) * P + 1):(g1 * P), ((g1 - 1) * P + 1):(g1 * P)]
  } else {
    Vg1 <- OmegaHat[((g1 - 1) * P + 1):(g1 * P), ((g1 - 1) * P + 1):(g1 * P)]
    Vg2 <- OmegaHat[((g2 - 1) * P + 1):(g2 * P), ((g2 - 1) * P + 1):(g2 * P)]
    Cov12 <- OmegaHat[((g1 - 1) * P + 1):(g1 * P), ((g2 - 1) * P + 1):(g2 * P)]
    Vdiff <- Vg1 + Vg2 - 2 * Cov12
  }
  
  # Statistics
  v_norm <- norm_vec(v)
  ZTv <- t(Z) %*% v
  ZTv_norm <- norm_vec(ZTv)
  dir_ZTv <- dir_vec(ZTv)
  Omega_ZTv_norm <- sqrt(t(ZTv) %*% solve(Vdiff) %*% ZTv)
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

#' Quadratic Form of ||z_it(ϕ) - ClusterMean_g(ϕ)||² for Canonical Form
#'
#' Represent
#' \deqn{||[z(ϕ)]_{it} - T^{-1} ∑_{t} ∑_{j} w_j^{(m-1)}(g) [z(ϕ)]_{jt}||²}
#' as a quadratic function in ϕ.
#'
#' @param rq List from \code{required_quantities()} (matrix-based)
#' @param i Integer; unit index
#' @param g Integer; cluster index
#' @param last_centroids K x P matrix of centroids from previous k-means step
#' @param weighted_deltahat K-length vector of weighted δₖ
#'
#' @return A list with \code{quad}, \code{linear}, \code{constant} terms
#' @keywords internal
#'
norm_phi_canonical_panel <- function(rq, i, g,
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
  D_deltahat <- deltahat_i - weighted_deltahat[g]
  D_Z <- Rfast::eachrow(Z_i, last_centroids[g, ], "-")
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
#' Computes \eqn{\sum_{j=1}^N w_j^{(m-1)}(g) \hat{\delta}_j} for each cluster g.
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
  for (g in seq_len(K)) {
    indices <- which(last_cl == g)
    weighted[g] <- if (length(indices) > 0) {
      sum(v1[indices]) / length(indices)
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
#' @export
compute_S <- function(rq, estimated_k_means) {
  requireNamespace(intervals)
  
  # Extract estimation path
  cluster_path <- do.call(rbind, estimated_k_means$clusters)
  centroid_path <- estimated_k_means$centers
  M <- nrow(cluster_path)
  
  # Panel dimensions and identifiers
  N <- rq$N
  K <- rq$K
  ids <- rq$ids
  Tobs <- rq$Tobs
  
  # Init containers
  all_interval_lists <- list()
  init_ids <- estimated_k_means$random_init_obs
  init_clustering <- cluster_path[1, ]
  
  # Initialization stage
  for (i in seq_len(N)) {
    id_i <- ids[i]
    g_star <- init_clustering[i]
    g_star_quad <- norm_sq_phi_panel(rq, id_i, init_ids[g_star])
    for (g in seq_along(init_ids)) {
      g_quad <- norm_sq_phi_panel(rq, id_i, init_ids[g])
      coeffs <- minus_quad_ineq(g_star_quad, g_quad)
      interval <- solve_one_ineq_complement(coeffs$quad, coeffs$linear, coeffs$constant)
      all_interval_lists[[(i - 1) * length(init_ids) + g]] <- interval
    }
  }
  
  # Refinement iterations
  index_offset <- length(all_interval_lists)
  count <- 1
  if (M > 1) {
    for (m in 2:M) {
      current_cl <- cluster_path[m, ]
      previous_cl <- cluster_path[m - 1, ]
      last_centroids <- centroid_path[[m]]
      wdelta <- weighted_deltahat(rq, previous_cl)
      
      for (i in seq_len(N)) {
        id_i <- ids[i]
        g_star <- current_cl[i]
        g_star_quad <- norm_phi_canonical_panel(rq, id_i, g_star, last_centroids, wdelta)
        
        for (g in seq_len(K)) {
          g_quad <- norm_phi_canonical_panel(rq, id_i, g, last_centroids, wdelta)
          coeffs <- minus_quad_ineq(g_star_quad, g_quad)
          interval <- solve_one_ineq_complement(coeffs$quad, coeffs$linear, coeffs$constant)
          all_interval_lists[[index_offset + count]] <- interval
          count <- count + 1
        }
      }
    }
  }
  
  # Final intersection: compute S
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

# ----- functions for computing tail probabilities of truncated chi-squared distributions -----
# ----- full credit to Shuxiao Chen, the writer of the outference package  -----
#' A helper function for approximating normal tail probabilities
#'
#' For \eqn{Z ~ N(0, 1)}, we have the approximation
#'     \eqn{P(Z \ge z) \approx }\code{magicfun(z)*exp(-z^2/2)}.
#'
#' @keywords internal
#'
#' @param z, the number where the function is evaluated.
#'
#' @return This function returns the value of the function evaluated at \code{z}.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
magicfun = function(z){
  z2 <- z*z
  z3 <- z*z*z
  temp <- (z2 + 5.575192695 * z + 12.77436324) /
    (sqrt(2*pi) * z3 + 14.38718147*z2 + 31.53531977*z + 2*12.77436324)
  return(temp)
}

#' Make endpoints of intervals finite
#'
#' This function modifies a union of intervals with positive but possibly infinite endpoints
#'    into a union of intervals with positive and \emph{finite} endpoints, while ensuring
#'    the probability of a \eqn{N(0, 1)} falling into it numerically the same.
#'
#' @keywords internal
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of intervals with \emph{positive} but possibly infinite endpoints.
#'
#' @return This function returns an "Intervals" object or a matrix depending on the input.
finiteE <- function(E) {
  ind.inf <- which(E == Inf)
  if (length(ind.inf) == 0) return(E)
  # we know there are some infinite entries
  E.max <- max(E[-ind.inf])
  E[which(E == Inf)] <- max(10000, E.max * 2)
  return(E)
}

#' Make endpoints of intervals positive
#'
#' This function modifies a union of intervals with possibly negative enpoints
#'     into a union of intervals with \emph{positive} endpoints, while ensuring
#'    the probability of a \eqn{N(0, 1)} falling into it numerically the same.
#'
#' @keywords internal
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of intervals with \emph{positive} but possibly infinite endpoints.
#'
#' @return This function returns an "Intervals" object or a matrix depending on the input.
sortE <- function(E) {
  E.sorted <- lapply(1:nrow(E), function(i){
    temp <- as.numeric(E[i, ])
    if (temp[1] <= 0 & temp[2] <= 0) {
      return(sort(-temp))
    }
    if (temp[1] >= 0 & temp[2] >= 0) {
      return(sort(temp))
    }
    # we know temp[1] < 0, temp[2] > 0 OR temp[1] > 0, temp[2] < 0
    temp <- abs(temp)
    return(rbind(c(0, temp[1]), c(0, temp[2])))
  })
  E.sorted <- do.call(rbind, E.sorted)
  # in order to use the approximation, we translate Inf to a large number
  return(finiteE(E.sorted))
}

#' Comparison between two intervals
#'
#' This functions returns \code{TRUE} if and only if two intervals are the same.
#'
#' @keywords internal
#'
#' @param int1,int2 "Intervals" objects.
#'
#' @return This function returns the desired logical result.
isSameIntervals <- function(int1, int2) {
  
  # first make int1, int2 to the default order
  int1 <- intervals::reduce(int1)
  int2 <- intervals::reduce(int2)
  
  if (nrow(int1) != nrow(int2)) return(FALSE)
  
  # int1 and int2 has the same number of intervals
  
  if (sum(int1 != int2) > 0) return(FALSE)
  
  # int1 and int2 has the same elements
  return(TRUE)
}

#' Approximation of the ratio of two normal probabilities
#'
#' This function returns an approximation of \eqn{P(Z \in E1)/P(Z \in E2)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param E1,E2 "Intervals" objects or matrices where rows represents
#'     a union of intervals with \emph{positive and finite} endpoints.
#' @param scale scaling parameter.
#'
#' @return This function returns the value of the approximation.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
TNRatioApprox <- function(E1, E2, scale = NULL) {
  
  if (is.null(scale)) {
    temp <- (c(E1, E2))^2
    scale.grid <- stats::quantile(temp, probs = seq(0, 1, 0.2))
    
    for(scale in scale.grid) {
      temp <- TNRatioApprox(E1, E2, scale = scale)
      if (!is.na(temp)) {
        return(temp)
      }
      # if temp is NaN, proceed to the next loop
    }
    
    # if all scale.grid does not work, then return NaN
    return(NaN)
  }
  num1 <- magicfun(E1[, 1]) * exp(-(E1[, 1]^2 - scale)/2)
  num2 <- magicfun(E1[, 2]) * exp(-(E1[, 2]^2 - scale)/2)
  denom1 <- magicfun(E2[, 1]) * exp(-(E2[, 1]^2 - scale)/2)
  denom2 <- magicfun(E2[, 2]) * exp(-(E2[, 2]^2 - scale)/2)
  res <- sum(num1-num2)/sum(denom1-denom2)
  return(res)
}

#' Probability of a standard normal in a single interval
#'
#' This function returns \eqn{P(lo \le Z \le up)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param lo,up quantiles.
#'
#' @return This function returns the desired probability.
TNProbEachInt <- function(lo, up) {
  if (up == Inf) {
    return(stats::pnorm(lo, 0, 1, lower.tail = FALSE))
  }
  # we know up < Inf, want P(lo <= X <= up).
  # we want small - small (big - big will mask small numbers),
  
  try1 <- stats::pnorm(lo, 0, 1, lower.tail = FALSE) - stats::pnorm(up, 0, 1, lower.tail = FALSE)
  if (try1 != 0) return(try1)
  
  try2 <- stats::pnorm(up, 0, 1, lower.tail = TRUE) - stats::pnorm(lo, 0, 1, lower.tail = TRUE)
  return(try2)
  
}

#' Probability of a standard normal in a union of intervals
#'
#' This function returns \eqn{P(Z \in E)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#'
#' @return This function returns the desired probability.
TNProb <- function(E) {
  # sum cdf over each disjoint interval of E
  res <- sum(sapply(1:nrow(E), function(v) {
    return(TNProbEachInt(E[v, 1], E[v, 2]))
  }))
  return(res)
}

#' Survival function of truncated normal distribution
#'
#' This function returns the upper tail probability of a truncated normal distribution
#'     at quantile \code{q}.
#'
#' Let \eqn{X} be a normal random variable with \code{mean} and \code{sd}. Truncating
#'     \eqn{X} to the set \eqn{E} is equivalent to conditioning on \eqn{{X \in E}}. So this function
#'     returns \eqn{P(X \ge q | X \in E)}.
#'
#' @keywords internal
#'
#' @param q the quantile.
#' @param mean the mean parameter
#' @param sd the standard deviation
#' @param E the truncation set, an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#' @param approx should the approximation algorithm be used? Default is \code{FALSE},
#'     where the approximation is not used in the first place. But when the result is wacky,
#'     the approximation will be used.
#'
#' @return This function returns the value of the survival function evaluated at quantile \code{q}.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
TNSurv <- function(q, mean, sd, E, approx = FALSE) {
  # check if truncation is empty (i.e. truncated to the empty set)
  if (nrow(E) == 0) {
    stop("truncation set is empty")
  }
  
  # check if truncation is the whole real line
  if (isSameIntervals(E, intervals::Intervals(c(-Inf, Inf)))) {
    return(stats::pnorm(q, mean, sd, lower.tail = FALSE))
  }
  
  # E is not empty and is not the whole real line,
  # i.e. 0 < P(X in E) < 1
  
  # we want P(X > q | X in E) = P(X >= q AND X in E) / P(X in E)
  # {X >= q} = {Z >= (q-mean)/sd}
  # {X in E} = {Z in (E-mean)/sd}
  # Z ~ N(0, 1)
  q <- (q-mean)/sd
  E <- (E-mean)/sd
  mean <- 0
  sd <- 1
  q2 <- q*q
  region <- suppressWarnings(intervals::interval_intersection(E, intervals::Intervals(c(q, Inf))))
  # check if the result is 0 or 1
  if(nrow(region) == 0) return(0)
  if (isSameIntervals(E, region)) return(1)
  
  # transform region and E so that intervals have positive endpoints
  region <- sortE(region)
  E <- sortE(E)
  
  # we want P(Z in region) / P(Z in E)
  # try approximate calculation
  if (approx) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }
  
  # try exact calculation
  denom <- TNProb(E)
  num <- TNProb(region)
  
  if (denom < 1e-100 || num < 1e-100) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }
  
  # we know denom and num are both reasonably > 0
  
  res <- num / denom
  # force the result to lie in [0, 1]
  return(max(0, min(1, res)))
}

#' Approximation of the ratio of two chi-squared probabilities
#'
#' This function returns an approximation of \eqn{P(X \in E1)/P(X \in E2)}, where
#'     \eqn{X} is a central chi-squared random variable with \code{df} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df degree of freedom of the chi-squared random variable.
#' @param E1,E2 "Intervals" objects or matrices where rows represents
#'     a union of intervals with \emph{positive and finite} endpoints.
#'
#' @return This function returns the value of the approximation.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
#' @references Canal, Luisa. "A normal approximation for the chi-square distribution."
#'     Computational statistics & data analysis 48.4 (2005): 803-808.
TChisqRatioApprox <- function(df, E1, E2) {
  
  # the transform that makes x into a N(0, 1) r.v. such that
  # P(X >= x) = P(Z >= Chisq2N(x)), X ~ chisq(df), Z ~ N(0, 1)
  # this function can take either scaler, vector or matrix
  Chisq2N <- function(x, df, tol = 1e-6) {
    
    if (is.numeric(x) && length(x) == 1) {
      if (x <= tol) { # x <= 0
        return(-Inf)
      }
      if (x == Inf) {
        return(Inf)
      }
      # we know x > 0 and x is finite
      x <- (x/df)^(1/6) - (1/2) * (x/df)^(1/3) + (1/3) * (x/df)^(1/2)
      mu <- 5/6 - 1/(9*df) - 7/(648*df^2) + 25/(2187*df^3)
      sig <- sqrt(1/(18*df) + 1/(162*df^2) - 37/(11664*df^3))
      return((x-mu)/sig)
    }
    
    if (is.vector(x)) {
      return(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df))
    }
    
    if (is.matrix(x)) {
      return(structure(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df), dim = dim(x)))
    }
    
    return(intervals::Intervals())
    
  }
  
  E1 <- Chisq2N(E1, df)
  E1 <- sortE(E1) # notice that Chisq2N can be negative
  E2 <- Chisq2N(E2, df)
  E2 <- sortE(E2)
  
  # now we want P(Z in E1) / P(Z in E2), Z ~ N(0, 1)
  return(TNRatioApprox(E1, E2))
}