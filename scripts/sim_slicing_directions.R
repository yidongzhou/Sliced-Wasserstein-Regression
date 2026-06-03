## Simulation to test different numbers of slicing directions for SWW
## Tests L = 50, 100, 200 * d for d = 2, 5, with n = 200
## Focuses on Gaussian setting to answer reviewer question about choice of slicing directions

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
L_multipliers <- c(50, 100, 200)  # Will be multiplied by d
n <- 200  # Fixed sample size as requested
max_cores <- 60
mc_cores <- max(1, min(detectCores(), max_cores))

# Model components (Gaussian setting)
alpha_fun <- function(x, d) {
  rep(x, d)
}

D_fun <- function(x, d) {
  diag(rep(x + 1, d))
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

run_one_mc <- function(run_id, d, L, n, xOut, mu_true, Sigma_true) {
  # Set seed for reproducibility within each run
  set.seed(run_id)
  
  # Simulate predictors
  x <- runif(n, min = -0.5, max = 0.5)
  
  # Simulate distributions (Gaussian setting only)
  y <- vector("list", n)
  N <- d * 100  # Sample size per distribution
  
  for (i in seq_len(n)) {
    xi <- x[i]
    mu <- alpha_fun(xi, d) + rnorm(d)
    Sigma <- rWishart(1, df = d + 1, Sigma = D_fun(xi, d))[,,1]
    
    # Generate Gaussian data
    y[[i]] <- mvrnorm(n = N, mu = mu, Sigma = Sigma)
  }
  
  # Run SWWR regression with specified number of directions
  res_swwr <- swwr(
    y = y,
    x = x,
    xOut = xOut,
    optns = list(method = "global", verbose = FALSE, cores = 1, L = L)
  )
  
  # Compute sliced Wasserstein errors for SWWR
  err_swwr <- vapply(seq_len(nOut), function(j) {
    Qpred_j <- t(sapply(res_swwr$sliced_quantiles, function(Q) Q[j, ]))  # L by M
    sliced_wasserstein_distance(mu = mu_true[j, ], Sigma = Sigma_true[[j]],
                                Qpred = Qpred_j,
                                directions = res_swwr$optns$directions)
  }, numeric(1))
  
  # Progress (printed from worker; order may be non-monotone)
  if (run_id %% 10 == 0) {
    cat(sprintf("Run %d completed (d=%d, L=%d)\n", run_id, d, L))
  }
  
  # Return average squared error across nOut samples
  mean(err_swwr^2)
}

# Main simulation loop
nOut <- 101
xOut <- seq(-0.5, 0.5, length.out = nOut)

for (d in dims) {
  cat(sprintf("\n=== Starting simulations for d = %d ===\n", d))
  
  # True regression parameters: E[mu|X=x] = alpha(x);
  # E[Sigma|X=x] = (d+1) * D(x)
  mu_true <- t(sapply(xOut, function(xi) alpha_fun(xi, d)))
  Sigma_true <- lapply(seq_len(nOut), function(i) (d + 1) * D_fun(xOut[i], d))
  
  # Collect results in a num_runs x (#L) matrix for this d
  run_errors_mat <- matrix(NA_real_, nrow = num_runs, ncol = length(L_multipliers))
  colnames(run_errors_mat) <- paste0("L=", (L_multipliers * d))
  
  # Store runtime for each L
  runtime_vec <- rep(NA_real_, length(L_multipliers))
  names(runtime_vec) <- paste0("L=", (L_multipliers * d))
  
  for (j in seq_along(L_multipliers)) {
    L <- L_multipliers[j] * d
    cat(sprintf("\n--- Testing L = %d (d = %d) ---\n", L, d))
    cat(sprintf("Starting %d runs with %d cores...\n", num_runs, mc_cores))
    
    # Time the simulation for this L
    start_time <- Sys.time()
    run_errors <- mcmapply(
      FUN = function(rid) run_one_mc(rid, d = d, L = L, n = n, 
                                     xOut = xOut, mu_true = mu_true, Sigma_true = Sigma_true),
      rid = seq_len(num_runs),
      mc.cores = mc_cores,
      mc.preschedule = FALSE
    )
    end_time <- Sys.time()
    runtime_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    # Store as a column in the matrix
    run_errors_mat[, j] <- as.numeric(run_errors)
    runtime_vec[j] <- runtime_seconds
    cat(sprintf("Completed L = %d (d = %d) successfully in %.1f seconds\n", L, d, runtime_seconds))
  }
  
  # Save intermediate results
  out_file <- swr_path("data", sprintf("slicing_directions_d%d.rds", d))
  runtime_file <- swr_path("data", sprintf("slicing_directions_runtime_d%d.rds", d))
  saveRDS(run_errors_mat, out_file)
  saveRDS(runtime_vec, runtime_file)
  cat(sprintf("Saved intermediate results to %s and %s\n", out_file, runtime_file))
}

# overall summary across dimensions (if files exist)
cat("\n=== SUMMARY (compact) ===\n")
for (d in dims) {
  f <- swr_path("data", sprintf("slicing_directions_d%d.rds", d))
  runtime_f <- swr_path("data", sprintf("slicing_directions_runtime_d%d.rds", d))
  if (file.exists(f) && file.exists(runtime_f)) {
    mat <- readRDS(f)
    runtime <- readRDS(runtime_f)
    means <- colMeans(mat, na.rm = TRUE)
    sds <- apply(mat, 2, sd, na.rm = TRUE)
    line <- paste(sprintf("%s mean=%.4g sd=%.4g time=%.1fs", colnames(mat), means, sds, runtime / 100), collapse = " | ")
    cat(sprintf("d=%d: %s\n", d, line))
  }
}
