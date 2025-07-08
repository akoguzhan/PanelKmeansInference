rm(list=ls())

# Required packages
library(cluster)   # For adjusted Rand index
library(mclust)    # For adjustedRandIndex
library(Matrix)    # For sparse matrix tools if needed
library(here)
library(parallel)
library(devtools)

# Load your packages
devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference")
devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA")

# Simulation parameters
set.seed(1)
N <- 120
Tobs <- 20
P <- 2
K_true <- 3
Kmax <- 5
n_sim <- 2000

# DGP function
simulate_panel_data <- function(N, Tobs, P, K, separation = 0) {
  group_sizes <- rep(floor(N/K), K)
  group_sizes[K] <- N - sum(group_sizes[1:(K-1)])
  group_labels <- rep(1:K, times = group_sizes)
  
  if (separation==0) {
    centers <- matrix(0, nrow = K, ncol = P)
  } else if (separation==1) {
    centers <- 0.2*rbind(c(4,2),c(2,1),c(1,0))
  }
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

alpha <- 0.05

# Parallel setup
n_cores <- max(1, detectCores(logical = FALSE) - 1)
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {
  library(mclust)
  devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/PanelKmeansInference", quiet = TRUE)
  devtools::load_all("C:/Users/og7051ak/Nextcloud/Codes for Papers/Tests with Unknown Clusters/Tools/clusteredEPA", quiet = TRUE)
})

clusterExport(cl, varlist = c("simulate_panel_data", "panel_homogeneity_test", "N", "Tobs", "P", "K_true", "alpha"),
              envir = environment())

sim_fun <- function(sim) {
  # separation = 1 power
  # separation = 0 size
  data <- simulate_panel_data(N, Tobs, P, K_true, separation = 1)
  Z <- data$Z
  id <- data$id
  time <- data$time
  true_cl <- data$true_cluster
  
  hom_test <- panel_homogeneity_test(
    Z = Z,
    id = id,
    time = time,
    K = 3,
    pcombine_fun = "Genmean_neq",
    r = -20,
    n_cores = 1 # Each worker uses 1 core
  )
  hom_pval <- hom_test$pval
  c(
    hom_pval,
  )
}

results <- parLapply(cl, 1:n_sim, sim_fun)
stopCluster(cl)
results <- do.call(rbind, results)

# Summary
cat("Mean <0.05 (hom_pval): ", mean(results[,1]<0.05), "\n")
cat("Mean <0.05 (grid_iu): ", mean(results[,2]<0.05), "\n")
cat("Mean <0.05 (iu_test): ", mean(results[,3]<0.05), "\n")
cat("Mean <0.05 (cauchy_test): ", mean(results[,4]<0.05), "\n")
cat("Mean <0.05 (simes_test): ", mean(results[,5]<0.05), "\n")
cat("Mean <0.05 (pairwise): ", mean(results[,6:8]<0.05), "\n")