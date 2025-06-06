rm(list=ls())

# Required packages
library(cluster)   # For adjusted Rand index
library(mclust)    # For adjustedRandIndex
library(Matrix)    # For sparse matrix tools if needed
library(here)
library(parallel)

devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")

# Simulation parameters
set.seed(1234)
N <- 150
Tobs <- 40
P <- 2
K_true <- 3
Kmax <- 5
n_sim <- 200

# DGP function
simulate_panel_data <- function(N, Tobs, P, K, separation = 0) {
  group_sizes <- rep(floor(N/K), K)
  group_sizes[K] <- N - sum(group_sizes[1:(K-1)])
  group_labels <- rep(1:K, times = group_sizes)
  
  centers <- matrix(rnorm(K * P, mean = 0, sd = separation), nrow = K)
  Z <- matrix(NA, N * Tobs, P)
  id_vec <- rep(1:N, each = Tobs)
  time_vec <- rep(1:Tobs, N)
  
  for (i in 1:N) {
    g <- group_labels[i]
    Z_i <- matrix(rep(centers[g, ], each = Tobs), nrow = Tobs) +
      matrix(rnorm(Tobs * P, sd = 1), nrow = Tobs)
    Z[((i - 1) * Tobs + 1):(i * Tobs), ] <- Z_i
  }
  
  list(Z = Z, id = id_vec, time = time_vec, true_cluster = group_labels)
}

# Performance metrics
rand_recovery <- function(true, est) {
  ari <- mclust::adjustedRandIndex(true, est)
  recovery <- sum(true == est) / length(true)
  return(c(rand = ari, recovery = recovery))
}

# Monte Carlo loop
results <- matrix(NA, nrow = n_sim, ncol = 5)
alpha <- 0.05

for (sim in 1:n_sim) {
  cat("Simulation", sim, "\n")
  data <- simulate_panel_data(N, Tobs, P, K_true)
  Z <- data$Z
  id <- data$id
  time <- data$time
  true_cl <- data$true_cluster
  
  # Selective inference homogeneity test
  hom_test <- panel_homogeneity_test(
    Z = Z,
    id = id,
    time = time,
    K = 3,
    pcombine_fun = "pmean",
    method = "A",
    r = 0,
    n_cores = (detectCores(logical = FALSE) - 1)
  )
  hom_pval <- hom_test$pvalue_combination
  results[sim,1] <- hom_pval
  results[sim,2] <- bonferroni_p(hom_test$pairwise_pvalues)
  results[sim,3:5] <- hom_test$pairwise_pvalues
}

# Summary
mean(results[,1]<0.05)
mean(results[,2]<0.05)
mean(results[,3:5]<0.05)