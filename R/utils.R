#' Fast Squared Euclidean Distance Matrix
#'
#' Computes the squared Euclidean distance between each row of matrix \code{A} and each row of matrix \code{B}.
#' This is a faster alternative to \code{proxy::dist()} for standard K-means clustering.
#'
#' @param A A numeric matrix of size \code{n1 x d}
#' @param B A numeric matrix of size \code{n2 x d}
#'
#' @return A numeric \code{n1 x n2} matrix where entry \code{(i,j)} is the squared Euclidean distance between row \code{i} of \code{A} and row \code{j} of \code{B}.
#' @examples
#' A <- matrix(rnorm(10 * 3), nrow = 10)
#' B <- matrix(rnorm(20 * 3), nrow = 20)
#' D <- squared_distance(A, B)
#' @keywords internal
squared_distance <- function(A, B) {
  A_sq <- rowSums(A^2)
  B_sq <- rowSums(B^2)
  outer(A_sq, B_sq, "+") - 2 * (A %*% t(B))
}

#' Compute Group Means by Cluster Membership
#'
#' @description Computes column means of a numeric matrix grouped by a vector of categorical indicators.
#'
#' @param X A numeric matrix of size N x P.
#' @param v A vector of group membership indicators of length N.
#'
#' @return A K x P matrix of group means, where K is the number of unique values in \code{v}.
#' @keywords internal
aggregate_matrix <- function(X, v) {
  if (!is.matrix(X)) stop("X must be a matrix.")
  v <- as.factor(v)
  counts <- as.numeric(table(v))
  group_sums <- rowsum(X, v)
  Xbar <- sweep(group_sums, 1, counts, FUN = "/")
  return(Xbar)
}

#' L2 Norm of a Vector
#'
#' @description Computes the Euclidean (L2) norm of a numeric vector.
#'
#' @param x A numeric vector.
#'
#' @return A numeric scalar representing the L2 norm of \code{x}.
#' @keywords internal
#' 
norm_vec = function(x) {
  l2_norm = sqrt(sum(x^2))
  return(l2_norm)
}

#' L2 Direction of a Vector
#'
#' @description Computes the normalized (L2 direction) vector.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector with L2 norm equal to 1.
#' @keywords internal
#' 
dir_vec = function(x) {
  l2_dir = x/norm_vec(x)
  return(l2_dir)
}

#' Check for Cluster Label Equivalence up to Permutation
#'
#' @description Compares two clustering vectors to determine if they differ only by a permutation.
#'
#' @param cl1 A vector of cluster labels.
#' @param cl2 Another vector of cluster labels.
#' @param K Number of clusters.
#'
#' @return Logical; \code{TRUE} if clusterings are identical up to permutation.
#' @keywords internal
#' 
same_cl = function(cl1,cl2,K) {
  tab = table(cl1,cl2)
  same_up_to_perm = sum(tab != 0) == K
  return(same_up_to_perm)
}

#' Convert Balanced Panel Data to Matrix Forms for EPA Testing and KMeans
#'
#' @param df A balanced panel data frame
#' @param id Character; name of unit identifier column
#' @param time Character; name of time identifier column
#' @param Z_names Character vector of variable names to extract
#'
#' @return A list with:
#' \describe{
#'   \item{Zbar}{T x P matrix of time-wise cross-sectional averages}
#'   \item{Z_panel}{NT x P matrix for clustering (unit-time stacked)}
#'   \item{id_numeric}{Numeric vector of unit identifiers (1 to N)}
#'   \item{time_numeric}{Numeric vector of time identifiers (1 to T)}
#' }
#' @keywords internal
#' 
panel_data_to_matrix <- function(df, id, time, Z_names) {
  df <- df[order(df[[id]], df[[time]]), ]
  
  # Convert id and time to numeric
  df$id_numeric <- as.numeric(factor(df[[id]]))
  df$time_numeric <- as.numeric(factor(df[[time]]))
  
  # T x P: time-averaged cross-sectional means
  Zbar_df <- aggregate(df[, Z_names], by = list(df[[time]]), FUN = mean)
  Zbar <- as.matrix(Zbar_df[, -1])
  
  # NT x P: stacked (unit-major) matrix
  Z_panel <- as.matrix(df[, Z_names])
  
  return(list(
    Zbar = Zbar,
    Z_panel = Z_panel,
    id_numeric = df$id_numeric,
    time_numeric = df$time_numeric
  ))
}
