## Monte Carlo simulation driver for Bayes-Hilbert baseline
## - Parallelizes outer MC; inner errors serial to avoid nesting
## - Saves per-run errors: Bayes global orig/trans, Bayes local orig/trans
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
dims <- c(2)  # Focus on d=2 only (bivariate)
n_vals <- c(50, 100, 200)
max_cores <- 60
mc_cores <- max(1, min(detectCores(), max_cores))

# Model components
# Global model functions
alpha_fun_global <- function(x, d) {
  rep(x, d)
}

D_fun_global <- function(x, d) {
  diag(rep(x + 1, d))
}

# Local model functions
alpha_fun_local <- function(x, d) {
  rep(sin(pi * x / 2) / 2, d)
}

D_fun_local <- function(x, d) {
  diag(rep(cos(pi * x / 2), d))
}

# Reproducible RNG across workers
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Sliced Wasserstein distance for density on grid (atoms + density matrix)
# Uses the more accurate approach: slice the 2D density to get 1D marginals, then compute quantiles
# Inputs:
#  - mu: d vector of true mean
#  - Sigma: d x d true covariance matrix
#  - atoms: K*L x d matrix of grid points
#  - density_mat: K x L matrix of predicted density values
#  - directions: L_dir x d matrix of direction vectors
#  - probs: optional probability grid for quantile comparison; default uses midpoints
# Output: scalar sliced W2 error for this observation
sliced_wasserstein_distance_bayes <- function(mu, Sigma, atoms, density_mat, directions, probs = NULL) {
  d <- length(mu)
  L_dir <- nrow(directions)
  M <- 200  # Number of quantiles to compute
  
  if (is.null(probs)) probs <- (seq_len(M) - 0.5) / M
  
  K <- nrow(density_mat)
  L <- ncol(density_mat)
  
  # Extract unique x and y coordinates from atoms to compute dx and dy
  # atoms are in expand.grid(x, y) order, so x varies fastest
  gx_unique <- unique(atoms[, 1])
  gy_unique <- unique(atoms[, 2])
  dx <- diff(gx_unique)[1]
  dy <- diff(gy_unique)[1]
  
  # Convert density matrix to vector matching atoms ordering
  density_vec <- as.vector(density_mat)  # column-major = expand.grid ordering
  
  # Compute quantiles for each direction by slicing the density
  Qpred <- matrix(0, L_dir, M)
  
  for (l in seq_len(L_dir)) {
    u <- directions[l, ]  # direction vector
    
    # Project all grid points along this direction: t = u' * (x, y)
    t_proj <- atoms %*% u  # (K*L) x 1
    
    # Get density values (already in correct order matching atoms)
    # Convert density to probability mass: pmf = density * dx * dy
    pmf_proj <- density_vec * dx * dy
    
    # Sort by projected values
    sorted_idx <- order(t_proj)
    t_sorted <- t_proj[sorted_idx]
    pmf_sorted <- pmf_proj[sorted_idx]
    
    # For a discrete grid, projections are typically distinct, so we can directly use sorted values
    # Normalize to ensure it's a proper PMF
    pmf_sorted <- pmf_sorted / sum(pmf_sorted)
    
    # Compute CDF from sorted projections
    cdf_marginal <- cumsum(pmf_sorted)
    
    # Compute quantiles by inverting the CDF
    Qpred[l, ] <- approx(cdf_marginal, t_sorted, probs, method = "linear",
                        yleft = min(t_sorted), yright = max(t_sorted))$y
  }
  
  # Vectorized means and variances across all directions
  mean_vec <- as.vector(directions %*% mu)                                     # L_dir
  var_vec <- rowSums((directions %*% Sigma) * directions)                      # L_dir
  sd_vec <- sqrt(pmax(var_vec, .Machine$double.eps))                           # L_dir
  
  # Build matrices for qnorm broadcast
  mean_mat <- matrix(mean_vec, nrow = L_dir, ncol = M)
  sd_mat <- matrix(sd_vec, nrow = L_dir, ncol = M)
  probs_mat <- matrix(probs, nrow = L_dir, ncol = M, byrow = TRUE)
  true_Q <- qnorm(probs_mat, mean = mean_mat, sd = sd_mat)                     # L_dir x M
  
  row_mse <- rowMeans((Qpred - true_Q)^2)                                      # L_dir
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
  
  # Generate directions for sliced Wasserstein distance
  directions <- runif_on_sphere(200, d)
  
  # Bayes baseline (global)
  gs <- 20  # K=20^2=400 for d=2
  res_bayes <- bayes_bivar_regression(y = y, x = x, xOut = xOut,
                                      optns = list(method = "global", grid_size = gs, verbose = FALSE))
  
  res_bayes_trans <- bayes_bivar_regression(y = yT, x = x, xOut = xOut,
                                            optns = list(method = "global", grid_size = gs, verbose = FALSE))
  
  # Compute errors
  nOut <- length(xOut)
  err_bayes <- numeric(nOut)
  err_bayes_trans <- numeric(nOut)
  
  for (j in seq_len(nOut)) {
    err_bayes[j] <- sliced_wasserstein_distance_bayes(
      mu = mu_true[j, ], Sigma = Sigma_true[[j]],
      atoms = res_bayes$atoms,
      density_mat = res_bayes$predicted_densities[[j]],
      directions = directions
    )
  }
  
  for (j in seq_len(nOut)) {
    err_bayes_trans[j] <- sliced_wasserstein_distance_bayes(
      mu = mu_true[j, ], Sigma = Sigma_true[[j]],
      atoms = res_bayes_trans$atoms,
      density_mat = res_bayes_trans$predicted_densities[[j]],
      directions = directions
    )
  }
  
  # Progress (printed from worker; order may be non-monotone)
  if (run_id %% 10 == 0) {
    cat(sprintf("%d runs completed\n", run_id))
  }
  
  # Return average error across nOut samples for both cases
  c(mean(err_bayes, na.rm = TRUE), mean(err_bayes_trans, na.rm = TRUE))
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
  
  # Generate directions for sliced Wasserstein distance
  directions <- runif_on_sphere(200, d)
  
  # Bayes baseline (local)
  gs <- 20  # K=20^2=400 for d=2
  res_bayes <- bayes_bivar_regression(y = y, x = x, xOut = xOut,
                                      optns = list(method = "local", grid_size = gs, verbose = FALSE, bw = bw_val))
  
  res_bayes_trans <- bayes_bivar_regression(y = yT, x = x, xOut = xOut,
                                            optns = list(method = "local", grid_size = gs, verbose = FALSE, bw = bw_val))
  
  # Compute errors
  nOut <- length(xOut)
  err_bayes <- numeric(nOut)
  err_bayes_trans <- numeric(nOut)
  
  for (j in seq_len(nOut)) {
    err_bayes[j] <- sliced_wasserstein_distance_bayes(
      mu = mu_true[j, ], Sigma = Sigma_true[[j]],
      atoms = res_bayes$atoms,
      density_mat = res_bayes$predicted_densities[[j]],
      directions = directions
    )
  }
  
  for (j in seq_len(nOut)) {
    err_bayes_trans[j] <- sliced_wasserstein_distance_bayes(
      mu = mu_true[j, ], Sigma = Sigma_true[[j]],
      atoms = res_bayes_trans$atoms,
      density_mat = res_bayes_trans$predicted_densities[[j]],
      directions = directions
    )
  }
  
  # Progress (printed from worker; order may be non-monotone)
  if (run_id %% 10 == 0) {
    cat(sprintf("%d runs completed\n", run_id))
  }
  
  # Return average error across nOut samples for both cases
  c(mean(err_bayes, na.rm = TRUE), mean(err_bayes_trans, na.rm = TRUE))
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
    
    # Use mclapply to get list, then convert to matrix
    run_errors_list <- mclapply(
      seq_len(num_runs),
      function(rid) run_one_mc_global(rid, d = d, N = N, n = n, 
                                      xOut = xOut, mu_true = mu_true, Sigma_true = Sigma_true),
      mc.cores = mc_cores,
      mc.preschedule = FALSE
    )
    # Convert list to matrix (num_runs x 2)
    run_errors <- do.call(rbind, run_errors_list)
    results_by_n[[i]] <- run_errors
    cat(sprintf("Completed global d=%d, n=%d successfully\n", d, n))
  }
  
  # Save results for this dimension
  out_file <- swr_path("data", sprintf("bayes_global_d%d.rds", d))
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
    
    # Use mclapply to get list, then convert to matrix
    run_errors_list <- mclapply(
      seq_len(num_runs),
      function(rid) run_one_mc_local(rid, d = d, N = N, n = n, 
                                     xOut = xOut, mu_true = mu_true, Sigma_true = Sigma_true),
      mc.cores = mc_cores,
      mc.preschedule = FALSE
    )
    # Convert list to matrix (num_runs x 2)
    run_errors <- do.call(rbind, run_errors_list)
    results_by_n[[i]] <- run_errors
    cat(sprintf("Completed local d=%d, n=%d successfully\n", d, n))
  }
  
  # Save results for this dimension
  out_file <- swr_path("data", sprintf("bayes_local_d%d.rds", d))
  saveRDS(results_by_n, out_file)
  cat(sprintf("Saved local results to %s\n", out_file))
}

# Global results
cat("\n--- GLOBAL RESULTS ---\n")
for (d in dims) {
  f <- swr_path("data", sprintf("bayes_global_d%d.rds", d))
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
        cat(sprintf("  %s (%d/%d valid) -> Bayes orig: %.4g(%.4g) | Bayes trans: %.4g(%.4g)\n",
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
  f <- swr_path("data", sprintf("bayes_local_d%d.rds", d))
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
        cat(sprintf("  %s (%d/%d valid) -> Bayes orig: %.4g(%.4g) | Bayes trans: %.4g(%.4g)\n",
                    nm, n_valid, nrow(error_mat), mu[1], sdv[1], mu[2], sdv[2]))
      }
    }
  } else {
    cat(sprintf("No local results file found for d=%d\n", d))
  }
}
