rm(list=ls())

# Required packages
library(cluster)   # For adjusted Rand index
library(mclust)    # For adjustedRandIndex
library(Matrix)    # For sparse matrix tools if needed
library(here)

source(here("Tools","PanelKmeansInference","estimation.R"))
source(here("Tools","PanelKmeansInference","utils.R"))

# Simulation parameters
set.seed(1234)
N <- 150
Tobs <- 40
P <- 2
K_true <- 3
Kmax <- 5
n_sim <- 100

# DGP function
simulate_panel_data <- function(N, Tobs, P, K, separation = 1) {
  group_sizes <- rep(floor(N/K), K)
  group_sizes[K] <- N - sum(group_sizes[1:(K-1)])
  group_labels <- rep(1:K, times = group_sizes)
  
  centers <- matrix(rnorm(K * P, mean = 0, sd = separation), nrow = K)
  if (K_true == 2) {
    centers <- rbind(c(4,2),c(2,1))
  } else {
    centers <- rbind(c(4,2),c(2,1),c(1,0),c(3,0))
  }
  
  Z <- matrix(NA, N * Tobs, P)
  id_vec <- rep(1:N, each = Tobs)
  
  for (i in 1:N) {
    g <- group_labels[i]
    Z_i <- matrix(rep(centers[g, ], each = Tobs), nrow = Tobs) +
      matrix(rnorm(Tobs * P, sd = 1), nrow = Tobs)
    Z[((i - 1) * Tobs + 1):(i * Tobs), ] <- Z_i
  }
  
  list(Z = Z, id = id_vec, true_cluster = group_labels)
}

# Performance metrics
rand_recovery <- function(true, est) {
  ari <- mclust::adjustedRandIndex(true, est)
  recovery <- sum(true == est) / length(true)
  return(c(rand = ari, recovery = recovery))
}

# Monte Carlo loop
results <- matrix(NA, nrow = n_sim, ncol = 4)
colnames(results) <- c("Rand", "Recovery", "BIC_Selects_Correct_K", "Khat")

for (sim in 1:n_sim) {
  print(sim)
  data <- simulate_panel_data(N, Tobs, P, K_true)
  Z <- data$Z
  id <- data$id
  true_cl <- data$true_cluster

  # Estimate with fixed K
  fit_fixed <- panel_kmeans_estimation(Z = Z, id = id, K = K_true, Ninit = 10, iter.max = 10)
  est_cl <- fit_fixed$final_cluster
  metrics <- rand_recovery(true_cl, est_cl)
  # Estimate with BIC
  fit_bic <- panel_kmeans_estimation(Z = Z, id = id, Kmax = Kmax, Ninit = 10, iter.max = 10)
  Khat <- length(unique(fit_bic$final_cluster))

  results[sim, ] <- c(metrics, Khat == K_true, Khat)
}

# Summary
summary_df <- data.frame(
  Rand_Index = mean(results[, "Rand"]),
  Recovery_Rate = mean(results[, "Recovery"]),
  BIC_Correct_Rate = mean(results[, "BIC_Selects_Correct_K"]),
  Avg_Khat = mean(results[, "Khat"])
)
print(summary_df)
