## Monte Carlo simulation driver for Fan & Müller (2024) baseline
## - Parallelizes outer MC; inner errors serial/subsampled to avoid nesting
## - Saves per-run errors: FM global orig/trans, FM local orig/trans
## - Focus on d=2 only for computational efficiency

library(parallel)
library(MASS)
library(T4transport)

.bootstrap_candidates <- c(file.path("scripts", "_bootstrap.R"), "_bootstrap.R")
.script_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[startsWith(commandArgs(trailingOnly = FALSE), "--file=")][1])
if (!is.na(.script_file)) {
  .bootstrap_candidates <- c(
    .bootstrap_candidates,
    file.path(dirname(normalizePath(.script_file, mustWork = FALSE)), "_bootstrap.R")
  )
}
.bootstrap_file <- .bootstrap_candidates[file.exists(.bootstrap_candidates)][1]
if (is.na(.bootstrap_file)) {
  stop("Cannot locate scripts/_bootstrap.R. Run from the project root or scripts directory.")
}
source(.bootstrap_file)
rm(.bootstrap_candidates, .bootstrap_file, .script_file)

# Configuration
num_runs <- 100
dims <- c(2)  # Focus on d=2 only for computational efficiency
n_vals <- c(50, 100, 200)
max_cores <- 60
mc_cores <- max(1, min(detectCores(), max_cores))

# Model components
# Global model functions (original from sim_fm.R)
alpha_fun_global <- function(x, d) {
  rep(x, d)
}

D_fun_global <- function(x, d) {
  diag(rep(x + 1, d))
}

# Local model functions (from sim_local.R)
alpha_fun_local <- function(x, d) {
  rep(sin(pi * x / 2) / 2, d)
}

D_fun_local <- function(x, d) {
  diag(rep(cos(pi * x / 2), d))
}

# Reproducible RNG across workers
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Sliced Wasserstein error for discrete PMF against true Gaussian regression
# Inputs:
#  - mu: d vector of true mean
#  - Sigma: d x d true covariance matrix
#  - atoms: K x d matrix of grid points
#  - pmf: K-length vector of probability masses
#  - directions: L x d matrix of direction vectors
#  - probs: optional probability grid for quantile comparison; default uses midpoints
# Output: scalar sliced W2 error for this observation
sliced_wasserstein_distance_fm <- function(mu, Sigma, atoms, pmf, directions, probs = NULL) {
  L <- nrow(directions)
  M <- 200  # Number of quantiles to compute
  
  if (is.null(probs)) probs <- (seq_len(M) - 0.5) / M
  
  # Project atoms along each direction
  projected_atoms <- atoms %*% t(directions)  # K x L
  
  # Compute quantiles for each direction
  Qpred <- matrix(0, L, M)
  for (l in seq_len(L)) {
    # Sort projected atoms and corresponding probabilities
    sorted_idx <- order(projected_atoms[, l])
    sorted_atoms <- projected_atoms[sorted_idx, l]
    sorted_pmf <- pmf[sorted_idx]
    
    # Compute cumulative distribution function
    cdf <- cumsum(sorted_pmf)
    
    # Interpolate to get quantiles
    Qpred[l, ] <- approx(cdf, sorted_atoms, probs, method = "linear", 
                        yleft = min(sorted_atoms), yright = max(sorted_atoms))$y
  }
  
  # Vectorized means and variances across all directions
  mean_vec <- as.vector(directions %*% mu)                                     # L
  var_vec <- rowSums((directions %*% Sigma) * directions)                      # L
  sd_vec <- sqrt(pmax(var_vec, .Machine$double.eps))                           # L
  
  # Build matrices for qnorm broadcast
  mean_mat <- matrix(mean_vec, nrow = L, ncol = M)
  sd_mat <- matrix(sd_vec, nrow = L, ncol = M)
  probs_mat <- matrix(probs, nrow = L, ncol = M, byrow = TRUE)
  true_Q <- qnorm(probs_mat, mean = mean_mat, sd = sd_mat)                     # L x M
  
  row_mse <- rowMeans((Qpred - true_Q)^2)                                      # L
  sqrt(mean(row_mse))
}


run_one_mc_global <- function(run_id, d, N, n, xOut, mu_true, Sigma_true) {
  # Set seed for reproducibility within each run
  set.seed(run_id)
  
  # Simulate predictors
  x <- runif(n, min = -0.5, max = 0.5)
  
  # Simulate distributions and apply transport map
  y <- vector("list", n)
  yT <- vector("list", n)
  
  for (i in seq_len(n)) {
    xi <- x[i]
    mu <- (alpha_fun_global(xi, d) + rnorm(d))
    S <- rWishart(1, df = d + 1, Sigma = D_fun_global(xi, d))[,,1]
    
    # Generate original data
    y[[i]] <- mvrnorm(n = N, mu = mu, Sigma = S)
    
    # Generate transport map parameter for this sample (same k for all coordinates)
    k <- sample(c(-2L, -1L, 1L, 2L), size = 1)
    
    # Apply transport map coordinate-wise with the same k
    yT[[i]] <- sweep(y[[i]], 2, rep(k, d), function(z, k) z - sin(k * z) / abs(k))
  }
  
  # Fan & Müller baseline (global entropic W2 barycenters)
  gs <- 20  # K=20^2=400 for d=2
  res_fm <- conditional_wb(y = y, x = x, xOut = xOut, 
                           optns = list(method = "global", verbose = FALSE, grid_size = gs, cores = 1)
                           )
  res_fm_trans <- conditional_wb(y = yT, x = x, xOut = xOut, 
                                 optns = list(method = "global", verbose = FALSE, grid_size = gs, cores = 1)
                                 )
  
  # Generate directions for sliced Wasserstein distance
  directions <- runif_on_sphere(200, d)
  
  # Compute sliced Wasserstein errors for FM original data
  err_fm <- vapply(seq_len(length(xOut)), function(j) {
    sliced_wasserstein_distance_fm(mu = mu_true[j, ], Sigma = Sigma_true[[j]],
                                  atoms = res_fm$atoms, pmf = res_fm$pred_bary[j, ],
                                  directions = directions)
  }, numeric(1))
  
  # Compute sliced Wasserstein errors for FM transported data
  err_fm_trans <- vapply(seq_len(length(xOut)), function(j) {
    sliced_wasserstein_distance_fm(mu = mu_true[j, ], Sigma = Sigma_true[[j]],
                                  atoms = res_fm_trans$atoms, pmf = res_fm_trans$pred_bary[j, ],
                                  directions = directions)
  }, numeric(1))
  
  # Progress (printed from worker; order may be non-monotone)
  if (run_id %% 10 == 0) {
    cat(sprintf("%d runs completed\n", run_id))
  }
  
  # For this run, return average error across nOut samples for both cases
  c(mean(err_fm), mean(err_fm_trans))
}

run_one_mc_local <- function(run_id, d, N, n, xOut, mu_true, Sigma_true) {
  # Set seed for reproducibility within each run
  set.seed(run_id)
  
  # Simulate predictors
  x <- runif(n, min = -0.5, max = 0.5)
  
  # Simulate distributions and apply transport map
  y <- vector("list", n)
  yT <- vector("list", n)
  
  for (i in seq_len(n)) {
    xi <- x[i]
    mu <- (alpha_fun_local(xi, d) + rnorm(d))
    S <- rWishart(1, df = d + 1, Sigma = D_fun_local(xi, d))[,,1]
    
    # Generate original data
    y[[i]] <- mvrnorm(n = N, mu = mu, Sigma = S)
    
    # Generate transport map parameter for this sample (same k for all coordinates)
    k <- sample(c(-2L, -1L, 1L, 2L), size = 1)
    
    # Apply transport map coordinate-wise with the same k
    yT[[i]] <- sweep(y[[i]], 2, rep(k, d), function(z, k) z - sin(k * z) / abs(k))
  }
  
  # Bandwidth for local regression: n^{-0.2}/4 (same as sim_local.R)
  bw_val <- n^(-0.2) / 4
  
  # Fan & Müller baseline (local entropic W2 barycenters)
  gs <- 20  # K=20^2=400 for d=2
  res_fm <- conditional_wb(y = y, x = x, xOut = xOut, 
                           optns = list(method = "local", verbose = FALSE, grid_size = gs, cores = 1, bw = bw_val)
                           )
  res_fm_trans <- conditional_wb(y = yT, x = x, xOut = xOut, 
                                 optns = list(method = "local", verbose = FALSE, grid_size = gs, cores = 1, bw = bw_val)
                                 )
  
  # Generate directions for sliced Wasserstein distance
  directions <- runif_on_sphere(200, d)
  
  # Compute sliced Wasserstein errors for FM original data
  err_fm <- vapply(seq_len(length(xOut)), function(j) {
    sliced_wasserstein_distance_fm(mu = mu_true[j, ], Sigma = Sigma_true[[j]],
                                  atoms = res_fm$atoms, pmf = res_fm$pred_bary[j, ],
                                  directions = directions)
  }, numeric(1))
  
  # Compute sliced Wasserstein errors for FM transported data
  err_fm_trans <- vapply(seq_len(length(xOut)), function(j) {
    sliced_wasserstein_distance_fm(mu = mu_true[j, ], Sigma = Sigma_true[[j]],
                                  atoms = res_fm_trans$atoms, pmf = res_fm_trans$pred_bary[j, ],
                                  directions = directions)
  }, numeric(1))
  
  # Progress (printed from worker; order may be non-monotone)
  if (run_id %% 10 == 0) {
    cat(sprintf("%d runs completed\n", run_id))
  }
  
  # For this run, return average error across nOut samples for both cases
  c(mean(err_fm), mean(err_fm_trans))
}

# Main simulation loop
nOut <- 101
xOut <- seq(-0.5, 0.5, length.out = nOut)

# Global simulation
cat("\n=== GLOBAL SIMULATION ===\n")
for (d in dims) {
  results_by_n <- vector("list", length = length(n_vals))
  names(results_by_n) <- paste0("n = ", n_vals)
  # True regression parameters for global model: E[mu|X=x] = alpha_global(x);
  # E[Sigma|X=x] = (d+1) * D_global(x)
  mu_true <- t(sapply(xOut, function(xi) alpha_fun_global(xi, d)))
  Sigma_true <- lapply(seq_len(nOut), function(i) (d + 1) * D_fun_global(xOut[i], d))
  
  for (i in seq_along(n_vals)) {
    n <- n_vals[i]
    N <- d * 100
    
    cat(sprintf("\n=== Global Monte Carlo: d=%d, n=%d ===\n", d, n))
    cat(sprintf("Starting %d runs with %d cores...\n", num_runs, mc_cores))
    
    # Use mcmapply to get matrix directly (num_runs x 2)
    run_errors <- mcmapply(
      FUN = function(rid) run_one_mc_global(rid, d = d, N = N, n = n, 
                                            xOut = xOut, mu_true = mu_true, Sigma_true = Sigma_true),
      rid = seq_len(num_runs),
      mc.cores = mc_cores,
      mc.preschedule = FALSE
    )
    results_by_n[[i]] <- t(run_errors)  # transpose to get num_runs x 2
    cat(sprintf("Completed global d=%d, n=%d successfully\n", d, n))
  }
  
  # Save results for this dimension
  out_file <- swr_path("data", sprintf("fm_global_d%d.rds", d))
  saveRDS(results_by_n, out_file)
  cat(sprintf("Saved global results to %s\n", out_file))
}

# Local simulation
cat("\n=== LOCAL SIMULATION ===\n")
for (d in dims) {
  results_by_n <- vector("list", length = length(n_vals))
  names(results_by_n) <- paste0("n = ", n_vals)
  # True regression parameters for local model: E[mu|X=x] = alpha_local(x);
  # E[Sigma|X=x] = (d+1) * D_local(x)
  mu_true <- t(sapply(xOut, function(xi) alpha_fun_local(xi, d)))
  Sigma_true <- lapply(seq_len(nOut), function(i) (d + 1) * D_fun_local(xOut[i], d))
  
  for (i in seq_along(n_vals)) {
    n <- n_vals[i]
    N <- d * 100
    
    cat(sprintf("\n=== Local Monte Carlo: d=%d, n=%d ===\n", d, n))
    cat(sprintf("Starting %d runs with %d cores...\n", num_runs, mc_cores))
    
    # Use mcmapply to get matrix directly (num_runs x 2)
    run_errors <- mcmapply(
      FUN = function(rid) run_one_mc_local(rid, d = d, N = N, n = n, 
                                           xOut = xOut, mu_true = mu_true, Sigma_true = Sigma_true),
      rid = seq_len(num_runs),
      mc.cores = mc_cores,
      mc.preschedule = FALSE
    )
    results_by_n[[i]] <- t(run_errors)  # transpose to get num_runs x 2
    cat(sprintf("Completed local d=%d, n=%d successfully\n", d, n))
  }
  
  # Save results for this dimension
  out_file <- swr_path("data", sprintf("fm_local_d%d.rds", d))
  saveRDS(results_by_n, out_file)
  cat(sprintf("Saved local results to %s\n", out_file))
}

# Global results
cat("\n--- GLOBAL RESULTS ---\n")
for (d in dims) {
  f <- swr_path("data", sprintf("fm_global_d%d.rds", d))
  if (file.exists(f)) {
    res <- readRDS(f)
    cat(sprintf("\nGlobal Dimension d=%d:\n", d))
    for (nm in names(res)) {
      error_mat <- res[[nm]]  # num_runs x 2 matrix
      if (all(is.na(error_mat))) {
        cat(sprintf("  %s: FAILED (all runs failed)\n", nm))
      } else {
        mu <- colMeans(error_mat, na.rm = TRUE)
        sdv <- apply(error_mat, 2, sd, na.rm = TRUE)
        n_valid <- sum(!is.na(error_mat[, 1]))
        cat(sprintf("  %s (%d/%d valid) -> FM orig: %.4g(%.4g) | FM trans: %.4g(%.4g)\n",
                    nm, n_valid, nrow(error_mat), mu[1], sdv[1], mu[2], sdv[2]))
      }
    }
  } else {
    cat(sprintf("No global results file found for d=%d\n", d))
  }
}

# Local results
cat("\n--- LOCAL RESULTS ---\n")
for (d in dims) {
  f <- swr_path("data", sprintf("fm_local_d%d.rds", d))
  if (file.exists(f)) {
    res <- readRDS(f)
    cat(sprintf("\nLocal Dimension d=%d:\n", d))
    for (nm in names(res)) {
      error_mat <- res[[nm]]  # num_runs x 2 matrix
      if (all(is.na(error_mat))) {
        cat(sprintf("  %s: FAILED (all runs failed)\n", nm))
      } else {
        mu <- colMeans(error_mat, na.rm = TRUE)
        sdv <- apply(error_mat, 2, sd, na.rm = TRUE)
        n_valid <- sum(!is.na(error_mat[, 1]))
        cat(sprintf("  %s (%d/%d valid) -> FM orig: %.4g(%.4g) | FM trans: %.4g(%.4g)\n",
                    nm, n_valid, nrow(error_mat), mu[1], sdv[1], mu[2], sdv[2]))
      }
    }
  } else {
    cat(sprintf("No local results file found for d=%d\n", d))
  }
}
