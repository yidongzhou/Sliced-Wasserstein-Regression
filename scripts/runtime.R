## Runtime comparison for SWW, SAW, FM, and Bayes methods
## - Compares global vs local approaches for SWWR, SAWR, FM, and Bayes
## - Tests sample sizes n = 50, 100, 200 with nOut = 101
## - FM and Bayes methods only tested for d=2 (too slow for d=5)
## - Reports timing for each method and configuration

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
dims <- c(2, 5)
n_vals <- c(50, 100, 200)
nOut <- 101
xOut <- seq(-0.5, 0.5, length.out = nOut)

# Model components (global variant)
alpha_fun_global <- function(x, d) {
  rep(x, d)
}

D_fun_global <- function(x, d) {
  diag(rep(x + 1, d))
}

# Model components (local variant)
alpha_fun_local <- function(x, d) {
  rep(sin(pi * x / 2) / 2, d)
}

D_fun_local <- function(x, d) {
  diag(rep(cos(pi * x / 2), d))
}

# Reproducible RNG
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Function to generate test data
generate_test_data <- function(n, d, N, alpha_fun, D_fun) {
  # Simulate predictors
  x <- runif(n, min = -0.5, max = 0.5)
  
  # Simulate distributions
  y <- vector("list", n)
  
  for (i in seq_len(n)) {
    xi <- x[i]
    mu <- (alpha_fun(xi, d) + rnorm(d))
    S <- rWishart(1, df = d + 1, Sigma = D_fun(xi, d))[,,1]
    
    # Generate data
    y[[i]] <- mvrnorm(n = N, mu = mu, Sigma = S)
  }
  
  return(list(x = x, y = y))
}

# Function to time a single method run
time_method <- function(method_name, method_type, data, xOut, d, n, N) {
  cat(sprintf("Timing %s (%s) - d=%d, n=%d, N=%d...\n", method_name, method_type, d, n, N))
  
  if (method_name == "SWWR") {
    if (method_type == "global") {
      optns <- list(method = "global", verbose = FALSE, cores = 10)
    } else {  # local
      bw_val <- n^(-0.2) / 4
      optns <- list(method = "local", verbose = FALSE, cores = 10, bw = bw_val)
    }
    
    start_time <- Sys.time()
    result <- swwr(
      y = data$y,
      x = data$x,
      xOut = xOut,
      optns = optns
    )
    end_time <- Sys.time()
    
  } else if (method_name == "SAWR") {
    if (method_type == "global") {
      optns <- list(method = "global", verbose = FALSE, cores = 10)
    } else {  # local
      bw_val <- n^(-0.2) / 4
      optns <- list(method = "local", verbose = FALSE, cores = 10, bw = bw_val)
    }
    
    start_time <- Sys.time()
    result <- sawr(
      y = data$y,
      x = data$x,
      xOut = xOut,
      optns = optns
    )
    end_time <- Sys.time()
    
  } else if (method_name == "FM") {
    if (method_type == "global") {
      optns <- list(method = "global", verbose = FALSE, cores = 10, grid_size = 20)
    } else {  # local
      bw_val <- n^(-0.2) / 4
      optns <- list(method = "local", verbose = FALSE, cores = 10, bw = bw_val, grid_size = 20)
    }
    
    start_time <- Sys.time()
    result <- conditional_wb(
      y = data$y,
      x = data$x,
      xOut = xOut,
      optns = optns
    )
    end_time <- Sys.time()
    
  } else if (method_name == "Bayes") {
    if (method_type == "global") {
      optns <- list(method = "global", verbose = FALSE, grid_size = 20)
    } else {  # local
      bw_val <- n^(-0.2) / 4
      optns <- list(method = "local", verbose = FALSE, bw = bw_val, grid_size = 20)
    }
    
    start_time <- Sys.time()
    result <- bayes_bivar_regression(
      y = data$y,
      x = data$x,
      xOut = xOut,
      optns = optns
    )
    end_time <- Sys.time()
  }
  
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf("  Completed in %.2f seconds\n", runtime))
  
  return(runtime)
}

# Main runtime comparison
cat("=== RUNTIME COMPARISON: SWW vs SAW vs FM vs Bayes ===\n")
cat("Testing global and local approaches for SWW, SAW, FM, and Bayes methods\n")
cat("Sample sizes: n =", paste(n_vals, collapse = ", "), "\n")
cat("Output points: nOut =", nOut, "\n")
cat("Note: FM and Bayes methods only tested for d=2 (too slow for d=5)\n\n")

# Initialize results storage
results <- list()

for (d in dims) {
  N <- d * 100  # Same as in simulations
  cat(sprintf("\n=== DIMENSION d = %d (N = %d) ===\n", d, N))
  
  results[[paste0("d", d)]] <- list()
  
  for (n in n_vals) {
    cat(sprintf("\n--- Sample size n = %d ---\n", n))
    
    # Generate test data for global and local models
    set.seed(123)  # Ensure reproducibility
    data_global <- generate_test_data(n, d, N, alpha_fun_global, D_fun_global)
    data_local <- generate_test_data(n, d, N, alpha_fun_local, D_fun_local)
    
    # Store results for this configuration
    config_results <- list()
    
    # Test SWWR Global
    config_results$swwr_global <- time_method("SWWR", "global", data_global, xOut, d, n, N)
    
    # Test SWWR Local
    config_results$swwr_local <- time_method("SWWR", "local", data_local, xOut, d, n, N)
    
    # Test SAWR Global
    config_results$sawr_global <- time_method("SAWR", "global", data_global, xOut, d, n, N)
    
    # Test SAWR Local
    config_results$sawr_local <- time_method("SAWR", "local", data_local, xOut, d, n, N)
    
    # Test FM and Bayes methods only for d=2 (too slow for d=5)
    if (d == 2) {
      # Test FM Global
      config_results$fm_global <- time_method("FM", "global", data_global, xOut, d, n, N)
      
      # Test FM Local
      config_results$fm_local <- time_method("FM", "local", data_local, xOut, d, n, N)
      
      # Test Bayes Global
      config_results$bayes_global <- time_method("Bayes", "global", data_global, xOut, d, n, N)
      
      # Test Bayes Local
      config_results$bayes_local <- time_method("Bayes", "local", data_local, xOut, d, n, N)
    }
    
    # Store results
    results[[paste0("d", d)]][[paste0("n", n)]] <- config_results
    
    # Print summary for this configuration
    cat(sprintf("\nSummary for d=%d, n=%d:\n", d, n))
    cat(sprintf("  SWWR Global: %.2f sec\n", config_results$swwr_global))
    cat(sprintf("  SWWR Local:  %.2f sec\n", config_results$swwr_local))
    cat(sprintf("  SAWR Global: %.2f sec\n", config_results$sawr_global))
    cat(sprintf("  SAWR Local:  %.2f sec\n", config_results$sawr_local))
    
    # Include FM and Bayes results only for d=2
    if (d == 2) {
      cat(sprintf("  FM Global:   %.2f sec\n", config_results$fm_global))
      cat(sprintf("  FM Local:    %.2f sec\n", config_results$fm_local))
      cat(sprintf("  Bayes Global: %.2f sec\n", config_results$bayes_global))
      cat(sprintf("  Bayes Local:  %.2f sec\n", config_results$bayes_local))
    }
    
    # Calculate speedup ratios
    swwr_ratio <- config_results$swwr_local / config_results$swwr_global
    sawr_ratio <- config_results$sawr_local / config_results$sawr_global
    swwr_vs_sawr_global <- config_results$sawr_global / config_results$swwr_global
    swwr_vs_sawr_local <- config_results$sawr_local / config_results$swwr_local
    
    cat(sprintf("  Local/Global ratio - SWWR: %.2fx, SAWR: %.2fx\n", swwr_ratio, sawr_ratio))
    cat(sprintf("  SAWR/SWWR ratio - Global: %.2fx, Local: %.2fx\n", swwr_vs_sawr_global, swwr_vs_sawr_local))
  }
}

# Save results
saveRDS(results, swr_path("data", "runtime_comparison.rds"))

# Print final summary table
cat("\n=== FINAL RUNTIME SUMMARY ===\n")
cat("Method\t\tType\t\td\tn\tTime (sec)\n")
cat("----------------------------------------\n")

for (d in dims) {
  for (n in n_vals) {
    config <- results[[paste0("d", d)]][[paste0("n", n)]]
    cat(sprintf("SWWR\t\tGlobal\t\t%d\t%d\t%.2f\n", d, n, config$swwr_global))
    cat(sprintf("SWWR\t\tLocal\t\t%d\t%d\t%.2f\n", d, n, config$swwr_local))
    cat(sprintf("SAWR\t\tGlobal\t\t%d\t%d\t%.2f\n", d, n, config$sawr_global))
    cat(sprintf("SAWR\t\tLocal\t\t%d\t%d\t%.2f\n", d, n, config$sawr_local))
    
    # Include FM and Bayes results only for d=2
    if (d == 2) {
      cat(sprintf("FM\t\tGlobal\t\t%d\t%d\t%.2f\n", d, n, config$fm_global))
      cat(sprintf("FM\t\tLocal\t\t%d\t%d\t%.2f\n", d, n, config$fm_local))
      cat(sprintf("Bayes\t\tGlobal\t\t%d\t%d\t%.2f\n", d, n, config$bayes_global))
      cat(sprintf("Bayes\t\tLocal\t\t%d\t%d\t%.2f\n", d, n, config$bayes_local))
    }
  }
}

cat("\n=== RUNTIME COMPARISON COMPLETE ===\n")
