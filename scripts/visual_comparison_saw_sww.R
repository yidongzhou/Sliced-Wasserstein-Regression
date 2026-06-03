### Visual Comparison of SAW and SWW Fitted 2D Distributions

library(ggplot2)
library(MASS)
library(ks)
library(parallel)
library(T4transport)
library(ggdensity)

# Source required functions
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
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
rm(.bootstrap_candidates, .bootstrap_file, .script_file)

# Sliced Wasserstein distance function (from sim_global.R)
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

# Set up simulation parameters
set.seed(123)
n <- 1000
d <- 2
N <- 1000  # Sample size per observation
xOut <- 0  # Single prediction point

# Model components (same as sim_global.R)
alpha_fun <- function(x, d) {
  rep(x, d)
}

D_fun <- function(x, d) {
  diag(rep(x + 1, d))
}

# Simulate training data
x <- runif(n, min = -0.5, max = 0.5)
y <- vector("list", n)

for (i in seq_len(n)) {
  xi <- x[i]
  mu <- (alpha_fun(xi, d) + rnorm(d))
  S <- rWishart(1, df = d + 1, Sigma = D_fun(xi, d))[,,1]
  y[[i]] <- mvrnorm(n = N, mu = mu, Sigma = S)
}

# True parameters at xOut = 0
mu_true <- alpha_fun(0, d)
Sigma_true <- (d + 1) * D_fun(0, d)

cat("=== Running SAW Regression ===\n")
# Run SAW regression
saw_result <- sawr(
  y = y,
  x = x,
  xOut = xOut,
  optns = list(
    method = "global", 
    density = TRUE  # Enable KDE for visualization
  )
)

cat("\n=== Running SWW Regression ===\n")
# Run SWW regression with reconstruction enabled
sww_result <- swwr(
  y = y,
  x = x,
  xOut = xOut,
  optns = list(
    method = "global", 
    reconstruct = TRUE,  # Enable density reconstruction
    limits = matrix(c(-3.5, 3.5, -3.5, 3.5), nrow = 2, ncol = 2)  # Set reconstruction domain for 95% mass
  )
)

# Create data frames for plotting
# True density: Generate true bivariate Gaussian density
true_x <- seq(-3.5, 3.5, length.out = 100)
true_y <- seq(-3.5, 3.5, length.out = 100)
true_grid <- expand.grid(x = true_x, y = true_y)
true_density <- mvtnorm::dmvnorm(as.matrix(true_grid), mean = mu_true, sigma = Sigma_true)
true_df <- data.frame(
  x = true_grid$x,
  y = true_grid$y,
  density = as.vector(true_density),
  method = "True"
)

# Extract predicted samples and densities
saw_samples <- saw_result$predicted_samples[[1]]  # M x 2 matrix
saw_density <- kde(saw_result$predicted_samples[[1]], xmin = rep(-3.5, 2), xmax = rep(3.5, 2), eval.points = true_grid)$estimate
saw_df <- data.frame(
  x = true_grid$x,
  y = true_grid$y,
  density = as.vector(saw_density),
  method = "SAW"
)

# Extract reconstructed density from SWW
sww_density <- sww_result$density[[1]]  # grid_size x grid_size matrix
sww_df <- data.frame(
  x = true_grid$x,
  y = true_grid$y,
  density = as.vector(sww_density),
  method = "SWW"
)

# Combine data for 1x3 comparison (order: SWW, True, SAW)
# Ensure proper factor ordering
sww_df$method <- factor(sww_df$method, levels = c("SWW", "True", "SAW"))
true_df$method <- factor(true_df$method, levels = c("SWW", "True", "SAW"))
saw_df$method <- factor(saw_df$method, levels = c("SWW", "True", "SAW"))
combined_df <- rbind(sww_df, true_df, saw_df)

# Create 1x3 density comparison plot
ggplot(combined_df, aes(x = x, y = y, z = density)) +
  geom_contour_filled() +
  facet_wrap(~ method, ncol = 3) +
  theme_void() +
  theme(
    text = element_text(size = 20),
    legend.position = "none"
  ) +
  coord_fixed(ratio = 1)
ggsave(swr_path("figures", "sww_true_saw_contour.pdf"), width = 12, height = 4)

# Create 1x3 contour plots
ggplot(combined_df, aes(x = x, y = y, z = density)) +
  geom_contour(color = "black", alpha = 0.7) +
  facet_wrap(~ method, ncol = 3) +
  theme_void() +
  theme(text = element_text(size = 20)) +
  coord_fixed(ratio = 1)

# Create 1x3 density heatmaps
ggplot(combined_df, aes(x = x, y = y, fill = density)) +
  geom_tile() +
  facet_wrap(~ method, ncol = 3) +
  scico::scale_fill_scico(limits = c(min(combined_df$density), max(combined_df$density)), palette = "lajolla", direction = -1) + 
  theme_void() +
  theme(
    text = element_text(size = 20),
    legend.position = "none"
  ) +
  coord_fixed(ratio = 1)
ggsave(swr_path("figures", "sww_true_saw_tile.pdf"), width = 12, height = 4)

# ===== QUANTITATIVE COMPARISON: SLICED WASSERSTEIN DISTANCE =====
cat("\n=== Computing Sliced Wasserstein Distances ===\n")

# Generate true samples for comparison
true_samples <- mvrnorm(n = N, mu = mu_true, Sigma = Sigma_true)

# Compute sliced Wasserstein distance for SWW
# Extract quantiles from SWW result (correctly format as L x M matrix)
sww_quantiles <- t(sapply(sww_result$sliced_quantiles, function(Q) Q[1, ]))  # L x M matrix
sww_directions <- sww_result$optns$directions      # L x d matrix

# Compute SWW sliced Wasserstein distance
sww_sw_distance <- sliced_wasserstein_distance(
  mu = mu_true, 
  Sigma = Sigma_true,
  Qpred = sww_quantiles,  # Already L x M
  directions = sww_directions
)

# Compute sliced Wasserstein distance for SAW
saw_sw_distance <- T4transport::swdist(
  saw_samples, 
  true_samples, 
  nproj = nrow(sww_directions)
)$distance

# Print quantitative comparison
cat(sprintf("Sliced Wasserstein Distance Results:\n"))
cat(sprintf("  SWW: %.6f\n", sww_sw_distance))
cat(sprintf("  SAW: %.6f\n", saw_sw_distance))
cat(sprintf("  Difference (SWW - SAW): %.6f\n", sww_sw_distance - saw_sw_distance))
cat(sprintf("  Relative difference: %.2f%%\n", 100 * (sww_sw_distance - saw_sw_distance) / saw_sw_distance))
