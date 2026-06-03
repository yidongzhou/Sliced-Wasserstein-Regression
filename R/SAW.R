### Sliced Wasserstein Regression for Multivariate Distributions (Direct)
library(parallel)
library(caret)
library(T4transport)
library(ks)

.swr_source_utils <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(frame) {
    file <- frame$ofile
    if (is.null(file)) NA_character_ else file
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles)]

  candidates <- c()
  if (length(ofiles) > 0) {
    candidates <- c(
      candidates,
      file.path(dirname(normalizePath(ofiles[length(ofiles)], mustWork = FALSE)), "utils.R")
    )
  }
  candidates <- c(candidates, file.path("R", "utils.R"), "utils.R")

  utils_file <- candidates[file.exists(candidates)][1]
  if (is.na(utils_file)) {
    stop("Cannot locate utils.R. Source R/load_swr.R or run from the project root.")
  }
  source(utils_file, chdir = TRUE)
}
.swr_source_utils()
rm(.swr_source_utils)

##### Helper Functions ############

# Inverse permutation utility
# Given a permutation vector perm, returns the inverse permutation inv_perm
# such that inv_perm[perm[i]] = i for all i
invPerm <- function(perm) {
  n <- length(perm)
  ip <- integer(n)
  ip[perm] <- seq_len(n)
  ip
}

# Compute sliced Wasserstein objective and gradient for a given candidate X
# This is the core computation for the SW barycenter optimization
# Args:
#   wgts: weight vector for each sample (normalized to sum to 1)
#   directions: L x d matrix of slicing directions (unit vectors)
#   X: M x d candidate barycenter matrix
#   df: list of n extended samples, each M x d
#   n: number of samples
#   cores: number of cores for parallelization
# Returns:
#   list with SW (scalar objective) and Deriv (M x d gradient matrix)
Step <- function(wgts, directions, X, df, n, cores = 1) {
  L <- nrow(directions)
  M <- nrow(X)
  d <- ncol(X)
  wgts <- wgts / sum(wgts)  # Normalize weights to sum to 1
  
  # Parallelize over directions if cores > 1
  if (cores > 1) {
    dir_results <- mclapply(seq_len(L), function(k) {
      u <- matrix(directions[k, ], ncol = 1)  # Current direction vector
      XT <- as.vector(X %*% u)               # Project X onto direction u
      orderXT <- order(XT)                   # Get sorting order for X
      inv_order <- invPerm(orderXT)          # Inverse permutation for alignment
      slice_SW <- 0                          # Accumulate objective for this direction
      slice_Deriv <- matrix(0, nrow = M, ncol = d)  # Accumulate gradient for this direction
      
      for (i in seq_along(df)) {
        yi <- df[[i]]                        # Current sample
        xiT <- as.vector(yi %*% u)          # Project sample onto direction u
        sorted_xiT <- sort(xiT)             # Sort projected sample
        xiT_matched <- sorted_xiT[inv_order] # Align to X's order using inverse permutation
        diff_vec <- XT - xiT_matched        # Difference vector
        Wass_slice <- sum(diff_vec^2)       # Wasserstein distance squared for this slice
        slice_SW <- slice_SW + wgts[i] * Wass_slice
        slice_Deriv <- slice_Deriv + wgts[i] * (matrix(diff_vec, ncol = 1) %*% t(u))
      }
      list(SW = slice_SW, Deriv = slice_Deriv)
    }, mc.cores = cores)
    
    # Aggregate results across all directions
    SW <- sum(sapply(dir_results, function(res) res$SW))
    Deriv <- Reduce("+", lapply(dir_results, function(res) res$Deriv))
  } else {
    # Serial fallback - same logic as parallel version
    SW <- 0
    Deriv <- matrix(0, nrow = M, ncol = d)
    for (k in seq_len(L)) {
      u <- matrix(directions[k, ], ncol = 1)
      XT <- X %*% u
      orderXT <- order(XT)
      inv_order <- invPerm(orderXT)
      for (i in seq_along(df)) {
        yi <- df[[i]]
        xiT <- yi %*% u
        sorted_xiT <- sort(xiT)
        xiT_matched <- sorted_xiT[inv_order]
        diff_vec <- as.vector(XT - xiT_matched)
        Wass_slice <- sum(diff_vec^2)
        SW <- SW + wgts[i] * Wass_slice
        Deriv <- Deriv + wgts[i] * (matrix(diff_vec, ncol = 1) %*% t(u))
      }
    }
  }
  
  list(SW = SW, Deriv = Deriv)
}

# Gradient descent optimizer for sliced Wasserstein barycenter
# Minimizes the weighted sliced Wasserstein objective using gradient descent
# Args:
#   df: list of n extended samples, each M x d
#   wgts: weight vector for each sample
#   directions: L x d matrix of slicing directions
#   X_init: M x d initial guess for barycenter
#   eta: step size for gradient descent
#   max_iter: maximum number of iterations
#   eps: convergence tolerance (relative change in objective)
#   cores: number of cores for parallelization
# Returns:
#   list with X (final barycenter) and conv (convergence flag)
SW_minimizer_finder <- function(df, wgts, directions, X_init, eta = 0.01, max_iter = 1000, eps = 0.001, cores = 1) {
  n <- length(df)
  X <- X_init
  sol <- Step(wgts, directions, X, df, n, cores)
  SW <- sol$SW
  Deriv <- sol$Deriv
  N_iters <- 0
  conv <- FALSE
  
  # Gradient descent loop
  while (N_iters < max_iter) {
    X_new <- X - eta * Deriv                    # Gradient descent step
    sol_new <- Step(wgts, directions, X_new, df, n, cores)
    SW_new <- sol_new$SW
    
    # Check for convergence (relative change in objective)
    if (abs(SW_new - SW) / max(SW, 1e-10) < eps) {
      X <- X_new
      SW <- SW_new
      Deriv <- sol_new$Deriv
      conv <- TRUE
      break
    }
    
    # Update for next iteration
    X <- X_new
    SW <- SW_new
    Deriv <- sol_new$Deriv
    N_iters <- N_iters + 1
  }
  
  list(X = X, conv = conv)
}

# Cross-validation bandwidth selection for local SAW regression
# Uses k-fold CV to select optimal bandwidth by minimizing SW distance on held-out data
# Args:
#   xin: n x p training predictor matrix
#   xOut: nOut x p prediction points
#   optns: options list containing kernel, eta, max_iter, eps, cores
#   y_extended: list of n extended samples for training
#   directions: L x d slicing directions
#   X_init: initial barycenter for warm start
# Returns:
#   optimal bandwidth vector of length p
bwCV_sawr <- function(xin, xOut, optns, y_extended, directions, X_init) {
  n <- nrow(xin)
  p <- ncol(xin)
  Kern <- kerFctn(optns$kernel)
  kFold <- ifelse(n > 30, 10, n)  # Use 10-fold CV or leave-one-out for small n
  folds <- createFolds(seq_len(n), kFold)
  
  # Objective function: mean CV error for given bandwidth vector
  objFctn <- function(bw_vec) {
    outCV <- numeric(length(folds))
    for (f in seq_along(folds)) {
      foldidx <- folds[[f]]
      testidx <- foldidx
      trainidx <- setdiff(seq_len(n), testidx)
      
      # Split data into train/test
      xin_train <- xin[trainidx, , drop = FALSE]
      xin_test <- xin[testidx, , drop = FALSE]
      df_train <- y_extended[trainidx]
      df_test <- y_extended[testidx]
      n_test <- length(testidx)
      
      # Compute kernel weights for test points using training data
      xMat_train <- t(apply(xin_test, 1, function(xt) {
        diffs <- sweep(xin_train, 2, xt, "-")  # Differences from training points
        u <- sweep(diffs, 2, bw_vec, "/")      # Normalize by bandwidth
        if (p == 1) {
          Kvec <- Kern(u[, 1])                 # 1D kernel evaluation
        } else {
          Kvec <- apply(u, 1, function(ui) prod(sapply(seq_len(p), function(dd) Kern(ui[dd]))))
        }
        return(Kvec)
      }))
      
      # Fit SW barycenters for test points using training data
      fit_SW <- lapply(seq_len(n_test), function(j) {
        wgts_j <- xMat_train[j, ]
        SW_minimizer_finder(df_train, wgts_j, directions, X_init, optns$eta, optns$max_iter, optns$eps, optns$cores)$X
      })
      
      # Compute SW distances between fitted and true test samples
      errors <- sapply(seq_len(n_test), function(j) T4transport::swdist(fit_SW[[j]], df_test[[j]], p = 2, nproj = nrow(directions))$distance)
      outCV[f] <- mean(errors)
    }
    mean(outCV)
  }
  
  # Set bandwidth search ranges using SetBwRange from utils.R
  if (is.null(optns$bwRange)) {
    bw_mins <- numeric(p)
    bw_maxs <- numeric(p)
    for (dd in seq_len(p)) {
      range_dd <- SetBwRange(xin[, dd], xOut[, dd], optns$kernel)
      bw_mins[dd] <- range_dd$min
      bw_maxs[dd] <- range_dd$max
    }
    optns$bwRange <- rbind(bw_mins, bw_maxs)
  }
  
  # Optimize bandwidth: 1D uses optimize(), multi-D uses optim()
  if (p == 1) {
    opt_res <- optimize(f = objFctn, interval = optns$bwRange[, 1])
    bw_opt <- opt_res$minimum
  } else {
    par_init <- colMeans(optns$bwRange)
    res_opt <- optim(par = par_init, fn = objFctn, lower = optns$bwRange[1, ], upper = optns$bwRange[2, ], method = "L-BFGS-B")
    bw_opt <- res_opt$par
  }
  bw_opt
}

##### Main Function ############

sawr <- function(y, x, xOut = NULL, optns = list()) {
  ## Sliced Wasserstein Regression for Multivariate Distributions (Direct)
  ##
  ## Inputs:
  ##   y: list of length n of empirical measures; each element is an n_i x d matrix (variable sample sizes allowed)
  ##   x: n x p matrix of predictors. If p = 1, a numeric vector of length n is also accepted.
  ##   xOut: optional nOut x p matrix of prediction points. If p = 1, a numeric vector of length nOut is also accepted.
  ##   optns: named list of options
  ##     - method: "global" (default) or "local" — Fréchet regression type
  ##     - directions: optional L x d matrix of slicing directions; each row must be non-zero. If provided,
  ##                   rows are normalized to unit length and L is set to nrow(directions).
  ##     - L: optional integer; number of directions to generate if 'directions' is not provided.
  ##          If both 'directions' and 'L' are missing, defaults to d * 100.
  ##     - cores: integer; number of parallel workers (default: max(1, detectCores()-1))
  ##     - verbose: logical; print progress (default TRUE)
  ##     - eta: step size for gradient descent (default 0.01)
  ##     - max_iter: maximum iterations for gradient descent (default 200)
  ##     - eps: convergence tolerance for gradient descent (default 0.01)
  ##     - kernel: kernel type for local regression: "gauss" (default), "rect", "epan", "quar", "gausvar"
  ##     - bw: bandwidth vector (length p) for local regression when method = "local"
  ##     - bwRange: 2 x p matrix of bandwidth ranges for selection (local method)
  ##     - density: logical; compute KDE if TRUE and d <= 5 (default FALSE)
  ##
  ## Outputs:
  ##   A list containing:
  ##   - predicted_samples: list of length nOut; each element is an M x d matrix of predicted empirical samples
  ##   - predicted_densities: optional list of length nOut; each element is a ks::kde object for the corresponding sample (if optns$density = TRUE and d <= 5)
  ##   - y, x, xOut, optns: the inputs (passed through). The field optns$directions contains the L x d
  ##     unit-norm direction matrix used.
  
  # Input validation
  if (is.null(y) || is.null(x)) {
    stop("Both y and x must be provided")
  }
  
  if (!is.list(y)) {
    stop("y must be a list of empirical measures")
  }
  
  # Allow x to be provided as a vector when p = 1
  if (!is.matrix(x)) {
    if (is.vector(x) || is.data.frame(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix, data.frame, or vector")
    }
  }
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  
  n <- length(y)
  p <- ncol(x)
  
  if (length(y) != nrow(x)) {
    stop("Number of observations in y and x must match")
  }
  
  # Determine dimension of distributions
  d <- ncol(y[[1]])
  
  # Handle xOut (mirror x handling when provided)
  if (is.null(xOut)) {
    xOut <- x
    nOut <- n
  } else {
    if (!is.matrix(xOut)) {
      if (is.vector(xOut) || is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else {
        stop("xOut must be a matrix, data.frame, or vector")
      }
    }
    if (is.null(dim(xOut))) {
      xOut <- matrix(xOut, ncol = 1)
    }
    nOut <- nrow(xOut)
  }
  
  # Defaults (tuned for speed)
  if (is.null(optns$method)) optns$method <- "global"
  if (is.null(optns$cores)) optns$cores <- max(1, parallel::detectCores() - 1)
  if (is.null(optns$verbose)) optns$verbose <- TRUE
  if (is.null(optns$eta)) optns$eta <- 0.01
  if (is.null(optns$max_iter)) optns$max_iter <- 200# 1000
  if (is.null(optns$eps)) optns$eps <- 0.01# 0.001
  if (is.null(optns$kernel)) optns$kernel <- "gauss"
  if (is.null(optns$density)) optns$density <- FALSE
  
  if (isTRUE(optns$verbose)) {
    cat("Starting Sliced Wasserstein Regression\n")
    cat("Dimensions: n =", n, ", p =", p, ", d =", d, "\n")
  }
  
  # Handle variable sample sizes by extending to common length
  sample_sizes <- sapply(y, nrow)
  
  if (length(unique(sample_sizes)) == 1) {
    # All sample sizes are exactly the same - no extension needed
    M <- sample_sizes[1]
    y_extended <- y
    if (isTRUE(optns$verbose)) {
      cat("All sample sizes equal:", M, "\n")
    }
  } else {
    # Sample sizes vary - need to extend
    size_range <- max(sample_sizes) - min(sample_sizes)
    size_ratio <- size_range / min(sample_sizes)
    
    if (size_ratio < 0.2) {
      # If sample sizes are close (within 20% of min), use min
      M <- min(sample_sizes)
      if (isTRUE(optns$verbose)) {
        cat("Sample sizes are close, using min length:", M, "\n")
      }
      y_extended <- lapply(y, function(yi) yi[sample(nrow(yi), M), , drop = FALSE])
    } else {
      # Otherwise use least common multiple (capped for practicality)
      M <- min(plcm(sample_sizes), n * max(sample_sizes), 5000)
      if (isTRUE(optns$verbose)) {
        cat("Sample sizes vary significantly, using least common multiple length:", M, "\n")
      }
      
      # Extend original data to common length M
      y_extended <- lapply(seq_along(y), function(i) {
        yi <- y[[i]]
        n_i <- nrow(yi)
        residual <- M %% n_i
        if (residual) {
          indices <- c(rep(1:n_i, each = M %/% n_i), sample(1:n_i, residual))
          yi[indices, , drop = FALSE]
        } else {
          rep_indices <- rep(1:n_i, each = M %/% n_i)
          yi[rep_indices, , drop = FALSE]
        }
      })
    }
  }
  
  if (isTRUE(optns$verbose)) {
    cat("Extended to common sample size M =", M, "\n")
  }
  
  # Step 1: Set or generate direction vectors
  if (!is.null(optns$directions)) {
    directions <- as.matrix(optns$directions)
    if (ncol(directions) != d) {
      stop("optns$directions must have d columns")
    }
    # normalize to unit length to ensure valid slicing directions
    norms <- sqrt(rowSums(directions * directions))
    if (any(norms == 0)) {
      stop("optns$directions contains zero rows; each direction must be non-zero for normalization")
    }
    directions <- directions / norms
    optns$L <- nrow(directions)
  } else {
    if (is.null(optns$L)) optns$L <- d * 100
    directions <- runif_on_sphere(optns$L, d)
  }
  
  if (isTRUE(optns$verbose)) {
    cat("Using", optns$L, "slicing directions on", d, "D unit sphere\n")
  }
  
  # Compute barycenter for warm start
  bary_wgts <- rep(1, n)
  X_init <- Reduce("+", y_extended) / length(y_extended)
  bary_res <- SW_minimizer_finder(y_extended, bary_wgts, directions, X_init, optns$eta, optns$max_iter, optns$eps, optns$cores)
  bary <- bary_res$X
  if (isTRUE(optns$verbose)) {
    cat("Barycenter computed\n")
  }
  
  # Bandwidth selection for local method
  if (optns$method == "local") {
    if (is.null(optns$bw)) {
      if (isTRUE(optns$verbose)) {
        cat("Selecting bandwidth via CV...\n")
      }
      optns$bw <- bwCV_sawr(x, xOut, optns, y_extended, directions, bary)
      if (isTRUE(optns$verbose)) {
        cat("Bandwidth selected:", paste(round(optns$bw, 4), collapse = " "), "\n")
      }
    }
  }
  
  # Compute weight matrix
  Kern <- kerFctn(optns$kernel)
  if (optns$method == "global") {
    xbar <- colMeans(x)
    Sigma <- cov(x) * (n - 1) / n
    invSigma <- solve(Sigma)
    xMat <- t(apply(xOut, 1, function(xo) {
      dxo <- xo - xbar
      dxi <- sweep(x, 1, xbar, "-")
      s <- as.vector(1 + dxi %*% (invSigma %*% dxo))
      return(s)
    }))
  } else {
    xMat <- t(apply(xOut, 1, function(xt) {
      diffs <- sweep(x, 2, xt, "-")
      u <- sweep(diffs, 2, optns$bw, "/")
      if (p == 1) {
        Kvec <- Kern(u[, 1])
      } else {
        Kvec <- apply(u, 1, function(ui) prod(sapply(seq_len(p), function(dd) Kern(ui[dd]))))
      }
      return(as.vector(Kvec))
    }))
  }
  
  if (isTRUE(optns$verbose)) {
    cat("Weights computed\n")
  }
  
  # Step 2: Perform weighted barycenter computations
  if (isTRUE(optns$verbose)) {
    cat("Performing weighted barycenter computations on", nOut, "points...\n")
  }
  
  outer_workers <- min(nOut, optns$cores)
  optns$inner_cores <- max(1, floor(optns$cores / outer_workers))
  if (optns$cores > 1) {
    fit <- mclapply(seq_len(nOut), function(i) {
      wgts_i <- xMat[i, ]
      SW_minimizer_finder(y_extended, wgts_i, directions, bary, optns$eta, optns$max_iter, optns$eps, optns$inner_cores)$X
    }, mc.cores = optns$cores)
  } else {
    fit <- lapply(seq_len(nOut), function(i) {
      wgts_i <- xMat[i, ]
      SW_minimizer_finder(y_extended, wgts_i, directions, bary, optns$eta, optns$max_iter, optns$eps, optns$inner_cores)$X
    })
  }
  
  if (isTRUE(optns$verbose)) {
    cat("Sliced Wasserstein Regression completed\n")
  }
  
  # Optional: Compute density estimates if requested and d <= 5
  predicted_densities <- NULL
  if (isTRUE(optns$density) && d <= 5) {
    if (isTRUE(optns$verbose)) {
      cat("Computing kernel density estimates...\n")
    }
    if (optns$cores > 1) {
      predicted_densities <- mclapply(fit, function(sample_mat) {
        ks::kde(sample_mat)
      }, mc.cores = optns$cores)
    } else {
      predicted_densities <- lapply(fit, function(sample_mat) {
        ks::kde(sample_mat)
      })
    }
    if (isTRUE(optns$verbose)) {
      cat("Density estimation completed\n")
    }
  } else if (isTRUE(optns$density) && d > 5) {
    warning("Density estimation skipped: dimension d > 5")
  }
  
  # Prepare output (store directions inside optns for a unified interface)
  optns$directions <- directions
  output <- list(
    predicted_samples = fit,
    y = y,
    x = x,
    xOut = xOut,
    optns = optns
  )
  if (!is.null(predicted_densities)) {
    output$predicted_densities <- predicted_densities
  }
  
  return(output)
}
