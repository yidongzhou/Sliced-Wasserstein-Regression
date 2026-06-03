## Monte Carlo simulation driver (local method for both SWWR and SAWR)
## - Parallelizes outer Monte Carlo runs across available cores
## - Sets swwr/sawr cores = 1 to avoid nested parallelism
## - Prints progress: "q runs are complete" from workers
## - Saves per-run sliced Wasserstein errors for each (d, N, n)
## - Tests both SWWR and SAWR methods using local approach

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
dims <- c(2, 5)
n_vals <- c(50, 100, 200)
max_cores <- 60
mc_cores <- max(1, min(detectCores(), max_cores))

# Model components (local variant - same as sim_swwr_local.R)
alpha_fun <- function(x, d) {
  rep(sin(pi * x / 2) / 2, d)
}

D_fun <- function(x, d) {
  diag(rep(cos(pi * x / 2), d))
}

# Reproducible RNG across workers
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Sliced Wasserstein error for sliced quantiles against true Gaussian regression
# Inputs:
#  - mu:      d vector of true mean
#  - Sigma:   d x d true covariance matrix
#  - Qpred: L x M matrix of predicted quantiles for this observation
#  - directions: L x d matrix of direction vectors
#  - probs: optional probability grid corresponding to columns of predicted quantiles; default uses midpoints
# Output: scalar sliced W2 error for this observation
sliced_wasserstein_distance <- function(mu, Sigma, Qpred, directions, probs = NULL) {
  d <- length(mu)
  L <- nrow(directions)
  M <- ncol(Qpred)
  if (is.null(probs)) probs <- (seq_len(M) - 0.5) / M
  
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

run_one_mc <- function(run_id, d, N, n, xOut, mu_true, Sigma_true) {
  # Set seed for reproducibility within each run
  set.seed(run_id)
  
  # Simulate predictors
  x <- runif(n, min = -0.5, max = 0.5)
  
  # Simulate distributions and apply transport map
  y <- vector("list", n)
  yT <- vector("list", n)
  
  for (i in seq_len(n)) {
    xi <- x[i]
    mu <- alpha_fun(xi, d) + rnorm(d)
    S <- rWishart(1, df = d + 1, Sigma = D_fun(xi, d))[,,1]
    
    # Generate original data
    y[[i]] <- mvrnorm(n = N, mu = mu, Sigma = S)
    
    # Generate transport map parameter for this sample (same k for all coordinates)
    k <- sample(c(-2L, -1L, 1L, 2L), size = 1)
    
    # Apply transport map coordinate-wise with the same k
    yT[[i]] <- sweep(y[[i]], 2, rep(k, d), function(z, k) z - sin(k * z) / abs(k))
  }
  
  # Bandwidth for local regression: n^{-0.2}/4 (same as sim_swwr_local.R)
  bw_val <- n^(-0.2) / 4
  
  # Run SWWR regression on original data (local)
  res_swwr <- swwr(
    y = y,
    x = x,
    xOut = xOut,
    optns = list(method = "local", verbose = FALSE, cores = 1, bw = bw_val)
  )
  
  # Run SWWR regression on transported data (local)
  res_swwr_trans <- swwr(
    y = yT,
    x = x,
    xOut = xOut,
    optns = list(method = "local", verbose = FALSE, cores = 1, bw = bw_val)
  )
  
  # Run SAWR regression on original data (local)
  res_sawr <- sawr(
    y = y,
    x = x,
    xOut = xOut,
    optns = list(method = "local", verbose = FALSE, cores = 1, bw = bw_val)
  )
  
  # Run SAWR regression on transported data (local)
  res_sawr_trans <- sawr(
    y = yT,
    x = x,
    xOut = xOut,
    optns = list(method = "local", verbose = FALSE, cores = 1, bw = bw_val)
  )
  
  # Compute sliced Wasserstein errors for SWWR original data
  err_swwr <- vapply(seq_len(nOut), function(j) {
    Qpred_j <- t(sapply(res_swwr$sliced_quantiles, function(Q) Q[j, ]))# L by M
    sliced_wasserstein_distance(mu = mu_true[j, ], Sigma = Sigma_true[[j]],
                                Qpred = Qpred_j,
                                directions = res_swwr$optns$directions)
  }, numeric(1))
  
  # Compute sliced Wasserstein errors for SWWR transported data
  err_swwr_trans <- vapply(seq_len(nOut), function(j) {
    Qpred_j <- t(sapply(res_swwr_trans$sliced_quantiles, function(Q) Q[j, ]))# L by M
    sliced_wasserstein_distance(mu = mu_true[j, ], Sigma = Sigma_true[[j]],
                                Qpred = Qpred_j,
                                directions = res_swwr_trans$optns$directions)
  }, numeric(1))
  
  # Generate true samples for SAWR comparison
  true_samples <- lapply(seq_len(nOut), function(i) {
    mvrnorm(n = N, mu = mu_true[i, ], Sigma = Sigma_true[[i]])
  })
  
  # Compute sliced Wasserstein errors for SAWR original data
  err_sawr <- vapply(seq_len(nOut), function(j) {
    T4transport::swdist(res_sawr$predicted_samples[[j]], true_samples[[j]], nproj = nrow(res_sawr$optns$directions))$distance
  }, numeric(1))
  
  # Compute sliced Wasserstein errors for SAWR transported data
  err_sawr_trans <- vapply(seq_len(nOut), function(j) {
    T4transport::swdist(res_sawr_trans$predicted_samples[[j]], true_samples[[j]], nproj = nrow(res_sawr_trans$optns$directions))$distance
  }, numeric(1))
  
  # Progress (printed from worker; order may be non-monotone)
  if (run_id %% 10 == 0) {
    cat(sprintf("%d runs completed\n", run_id))
  }
  
  # For this run, return average error across nOut samples for all four cases
  c(mean(err_swwr^2), mean(err_swwr_trans^2), mean(err_sawr^2), mean(err_sawr_trans^2))
}

# Main simulation loop
nOut <- 101
xOut <- seq(-0.5, 0.5, length.out = nOut)
for (d in dims) {
  results_by_n <- vector("list", length = length(n_vals))
  names(results_by_n) <- paste0("n = ", n_vals)
  # True regression parameters: E[mu|X=x] = alpha(x);
  # E[Sigma|X=x] = (d+1) * D(x)
  mu_true <- t(sapply(xOut, function(xi) alpha_fun(xi, d)))
  Sigma_true <- lapply(seq_len(nOut), function(i) (d + 1) * D_fun(xOut[i], d))
  
  for (i in seq_along(n_vals)) {
    n <- n_vals[i]
    N <- d * 100
    
    cat(sprintf("\n=== Monte Carlo: d=%d, n=%d ===\n", d, n))
    cat(sprintf("Starting %d runs with %d cores...\n", num_runs, mc_cores))
    
    # Use mcmapply to get matrix directly (num_runs x 4)
    run_errors <- mcmapply(
      FUN = function(rid) run_one_mc(rid, d = d, N = N, n = n, 
                                     xOut = xOut, mu_true = mu_true, Sigma_true = Sigma_true),
      rid = seq_len(num_runs),
      mc.cores = mc_cores,
      mc.preschedule = FALSE
    )
    results_by_n[[i]] <- t(run_errors)  # transpose to get num_runs x 4
    cat(sprintf("Completed d=%d, n=%d successfully\n", d, n))
  }
  
  # Save results for this dimension
  out_file <- swr_path("data", sprintf("swr_local_d%d.rds", d))
  saveRDS(results_by_n, out_file)
  cat(sprintf("Saved results to %s\n", out_file))
}

# Summarize saved outputs: mean and sd of errors for each case
cat("\n=== SIMULATION SUMMARY ===\n")
for (d in dims) {
  f <- swr_path("data", sprintf("swr_local_d%d.rds", d))
  if (file.exists(f)) {
    res <- readRDS(f)
    cat(sprintf("\nDimension d=%d:\n", d))
    for (nm in names(res)) {
      error_mat <- res[[nm]]  # num_runs x 4 matrix
      if (all(is.na(error_mat))) {
        cat(sprintf("  %s: FAILED (all runs failed)\n", nm))
      } else {
        mu <- colMeans(error_mat, na.rm = TRUE)
        sdv <- apply(error_mat, 2, sd, na.rm = TRUE)
        n_valid <- sum(!is.na(error_mat[, 1]))
        cat(sprintf("  %s (%d/%d valid) -> SWWR orig: mean=%.4g sd=%.4g | SWWR trans: mean=%.4g sd=%.4g | SAWR orig: mean=%.4g sd=%.4g | SAWR trans: mean=%.4g sd=%.4g\n",
                    nm, n_valid, nrow(error_mat), mu[1], sdv[1], mu[2], sdv[2], mu[3], sdv[3], mu[4], sdv[4]))
      }
    }
  } else {
    cat(sprintf("No results file found for d=%d\n", d))
  }
}
