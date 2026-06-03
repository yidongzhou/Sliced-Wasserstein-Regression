# conditional_wb.R: Implementation of Fan & Muller (2024) Conditional W2 Barycenters
# Custom Dykstra/Sinkhorn to handle negative weights (translation of MATLAB perform_dikstra_scaling.m)
# Dependencies: MASS (for dmvnorm, cov), caret (for createFolds)
library(MASS)
library(caret)
library(parallel)

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

#' Geometric Mean Proximal Operator for Barycenter Functional
#' 
#' Computes the geometric mean of probability distributions with potentially negative weights.
#' This handles the case where Fréchet regression weights can be negative.
#' 
#' @param p List of M probability vectors, each of length K
#' @param lambda Weight vector of length M (can contain negative values)
#' @return List of M probability vectors (all identical, representing the geometric mean)
geometric_mean <- function(p, lambda) {
  M <- length(p)
  K <- length(p[[1]])
  lambda <- lambda / sum(lambda)  # Normalize
  h <- rep(1, K)
  for (m in 1:M) {
    h <- h + lambda[m] * log(p[[m]] + 1e-12)  # Avoid log(0)
  }
  h <- exp(h)
  q <- replicate(M, h, simplify = FALSE)
  return(q)
}

#' Dykstra's Scaling Algorithm for Entropic Optimal Transport
#' 
#' Implements Dykstra's alternating projection algorithm for solving entropic optimal transport
#' problems. This is a custom implementation that handles negative weights in the barycenter problem.
#' 
#' @param K_func Function that applies the kernel matrix: K_func(vec) = kernel_matrix %*% vec
#' @param dims Vector of length 2: c(K, M) where K is number of atoms, M is number of measures
#' @param proxs List of two proximal operators: prox1 (geometric mean), prox2 (marginal constraint)
#' @param niter Maximum number of iterations (default: 500)
#' @param verb Logical; print progress if TRUE (default: FALSE)
#' @param tol Convergence tolerance (default: 1e-6)
#' @return List containing:
#'   - a, b: Scaling vectors
#'   - p_list: List of intermediate solutions
#'   - bary_r: Final barycenter probability vector
perform_dykstra_scaling <- function(K_func, dims, proxs, niter = 500, verb = FALSE, tol = 1e-6) {
  K <- dims[1]
  M <- dims[2]
  
  # Initialization
  U <- rep(1, K)
  a <- replicate(M, U, simplify = FALSE)
  b <- a
  u <- a; u1 <- a
  v <- a; v1 <- a
  
  p_list <- list(list(), list())
  
  for (i in 1:niter) {
    for (it in 1:2) {
      # Save previous
      a1 <- a; b1 <- b
      u2 <- u1; u1 <- u
      v2 <- v1; v1 <- v
      
      # Compute tilde a, tilde b, p1
      ta1 <- vector("list", M); tb1 <- vector("list", M)
      p1 <- vector("list", M)
      for (m in 1:M) {
        ta1[[m]] <- a1[[m]] * u2[[m]]
        tb1[[m]] <- b1[[m]] * v2[[m]]
        if (it == 1) {
          p1[[m]] <- ta1[[m]] * K_func(tb1[[m]])
        } else {
          p1[[m]] <- tb1[[m]] * K_func(ta1[[m]])
        }
      }
      
      # Apply prox
      p <- proxs[[it]](p1)
      
      # Update a, b, u, v
      for (m in 1:M) {
        if (it == 1) {
          # Avoid division by zero
          ratio <- p[[m]] / (p1[[m]] + 1e-12)
          a[[m]] <- ta1[[m]] * ratio
          b[[m]] <- tb1[[m]]
        } else {
          # Avoid division by zero
          ratio <- p[[m]] / (p1[[m]] + 1e-12)
          b[[m]] <- tb1[[m]] * ratio
          a[[m]] <- ta1[[m]]
        }
        # Avoid division by zero
        u[[m]] <- u2[[m]] * a1[[m]] / (a[[m]] + 1e-12)
        v[[m]] <- v2[[m]] * b1[[m]] / (b[[m]] + 1e-12)
      }
      
      # Save p1
      p_list[[it]] <- c(p_list[[it]], list(p1))
    }
  }
  
  # Extract final barycenter
  last_p1 <- p_list[[1]][[length(p_list[[1]])]]
  final_p <- proxs[[1]](last_p1)
  bary_r <- final_p[[1]]
  
  list(a = a, b = b, p_list = p_list, bary_r = bary_r)
}

#' Custom Entropic Barycenter with Negative Weight Support
#' 
#' Computes the entropic Wasserstein barycenter using Dykstra's algorithm, which can handle
#' negative weights that arise in Fréchet regression. This is an alternative to the standard
#' Sinkhorn algorithm that requires non-negative weights.
#' 
#' @param r_list List of n probability vectors (marginals), each of length K
#' @param weights Weight vector of length n (can contain negative values)
#' @param atoms K x d matrix of atom locations
#' @param lambda Entropic regularization parameter (default: 0.4)
#' @param max_iter Maximum number of Dykstra iterations (default: 500)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param verb Logical; print progress if TRUE (default: FALSE)
#' @return K-length vector of barycenter probabilities
entropic_barycenter <- function(r_list, weights, atoms, lambda, max_iter = 500, tol = 1e-6, verb = FALSE) {
  n <- length(r_list)  # M
  K <- length(r_list[[1]])
  
  # Precompute kernel matrix: exp(-C / lambda), C = ||atoms_i - atoms_j||^2
  atoms_sq <- rowSums(atoms^2)
  C <- outer(atoms_sq, rep(1, K)) + outer(rep(1, K), atoms_sq) - 2 * (atoms %*% t(atoms))
  kernel_matrix <- exp(-C / lambda)
  
  # Kernel function
  K_func <- function(vec) kernel_matrix %*% vec
  
  # Fixed marginals p0
  p0 <- r_list
  
  # Proxs
  prox1 <- function(p1) geometric_mean(p1, weights)
  prox2 <- function(p1) p0
  proxs <- list(prox1, prox2)
  
  # Run Dykstra
  result <- perform_dykstra_scaling(K_func, c(K, n), proxs, max_iter, verb, tol)
  
  # Return
  as.vector(result$bary_r)
}

#' Discretize Empirical Samples to Discrete Densities
#' 
#' Converts a list of empirical samples into discrete probability distributions on a regular grid.
#' This is necessary for computing entropic Wasserstein barycenters.
#' 
#' @param y List of n empirical measures; each element is an n_i x d matrix
#' @param grid_size Integer; grid resolution per dimension (default: 101)
#' @param bounds Optional 2 x d matrix specifying domain bounds; if NULL, uses 2.5%-97.5% quantiles
#' @return List containing:
#'   - atoms: K x d matrix of grid point locations
#'   - r_list: List of n probability vectors (discretized measures)
#'   - K: Total number of grid points (grid_size^d)
#'   - grid_size: Grid resolution per dimension
#'   - bounds: 2 x d matrix of domain bounds used
discretize_measures <- function(y, grid_size, bounds = NULL) {
  n <- length(y)
  d <- ncol(y[[1]])
  if (d == 0) stop("Response dimension d must be >=1")
  K <- grid_size^d
  if (K > 1e6) warning("K = ", K, " may be too large; reduce grid_size or d")
  
  # Data-driven bounds if not provided
  if (is.null(bounds)) {
    all_y <- do.call(rbind, y)
    bounds <- apply(all_y, 2, quantile, probs = c(0.025, 0.975))
    bounds <- apply(bounds, 2, range)
  }
  if (ncol(bounds) != d || nrow(bounds) != 2) stop("bounds must be 2 x d")
  
  grid_edges <- lapply(1:d, function(dim) seq(bounds[1, dim], bounds[2, dim], length.out = grid_size + 1))
  
  # Grid centers (atoms)
  g_indices <- arrayInd(1:K, rep(grid_size, d))
  atoms <- matrix(0, K, d)
  for (dim in 1:d) {
    atoms[, dim] <- (grid_edges[[dim]][g_indices[, dim]] + grid_edges[[dim]][g_indices[, dim] + 1]) / 2
  }
  
  r_list <- vector("list", n)
  for (i in 1:n) {
    yi <- y[[i]]
    n_i <- nrow(yi)
    if (n_i == 0) stop("Empty samples in y[[", i, "]]")
    
    # Multi-D binning
    bin_idx <- matrix(0, n_i, d)
    for (dim in 1:d) {
      bin_idx[, dim] <- findInterval(yi[, dim], grid_edges[[dim]])
      bin_idx[, dim] <- pmax(1, pmin(grid_size, bin_idx[, dim]))
    }
    
    # Linearize indices
    sub_d <- d - 1
    if (sub_d > 0) {
      sub_strides <- grid_size ^ ((d - 2):0)
      offset_sub <- apply(bin_idx[, 2:d, drop = FALSE], 1, function(row) sum(sub_strides * (row - 1)))
    } else {
      offset_sub <- rep(0, n_i)
    }
    flat_idx <- (bin_idx[, 1] - 1) * (grid_size ^ (d - 1)) + offset_sub + 1
    
    # Counts
    h <- tabulate(flat_idx, nbins = K)
    if (sum(h) == 0) h <- rep(1, K)
    r_list[[i]] <- h / sum(h)
  }
  
  list(atoms = atoms, r_list = r_list, K = K, grid_size = grid_size, bounds = bounds)
}

#' Cross-Validation Bandwidth Selection for Local FM Regression
#' 
#' Selects optimal bandwidth for local Fréchet regression using k-fold cross-validation.
#' The objective is to minimize the sliced Wasserstein distance between predicted and true distributions.
#' 
#' @param xin n x p training predictor matrix
#' @param xOut nOut x p prediction points
#' @param optns List of options containing kernel, lambda, max_iter, tol, grid_size, bounds
#' @param y List of n empirical measures for training
#' @return Optimal bandwidth vector of length p
bwCV_fm <- function(xin, xOut, optns, y) {
  n <- nrow(xin)
  p <- ncol(xin)
  Kern <- kerFctn(optns$kernel)
  kFold <- ifelse(n > 30, 10, n)
  folds <- createFolds(seq_len(n), kFold)
  
  # Objective function
  objFctn <- function(bw_vec) {
    outCV <- numeric(length(folds))
    for (f in seq_along(folds)) {
      foldidx <- folds[[f]]
      testidx <- foldidx
      trainidx <- setdiff(seq_len(n), testidx)
      
      # Split
      xin_train <- xin[trainidx, , drop = FALSE]
      xin_test <- xin[testidx, , drop = FALSE]
      y_train <- y[trainidx]
      y_test <- y[testidx]
      n_test <- length(testidx)
      
      # Discretize train
      disc_train <- discretize_measures(y_train, optns$grid_size, optns$bounds)
      atoms_train <- disc_train$atoms
      r_list_train <- disc_train$r_list
      K_train <- disc_train$K
      
      # Kernel weights
      xMat_train <- t(apply(xin_test, 1, function(xt) {
        diffs <- sweep(xin_train, 2, xt, "-")
        u <- sweep(diffs, 2, bw_vec, "/")
        if (p == 1) {
          Kvec <- Kern(u[, 1])
        } else {
          Kvec <- apply(u, 1, function(ui) prod(sapply(seq_len(p), function(dd) Kern(ui[dd]))))
        }
        return(Kvec)
      }))
      
      # Fit barycenters
      if (optns$cores > 1) {
        # Parallel computation for CV
        fit_fm <- mclapply(seq_len(n_test), function(j) {
          wgts_j <- xMat_train[j, ]
          if (sum(wgts_j) == 0) wgts_j <- rep(1/length(wgts_j), length(wgts_j))
          else wgts_j <- wgts_j / sum(wgts_j)
          entropic_barycenter(r_list_train, wgts_j, atoms_train, optns$lambda, optns$max_iter, optns$tol)
        }, mc.cores = min(optns$cores, n_test))
      } else {
        # Serial computation for CV
        fit_fm <- lapply(seq_len(n_test), function(j) {
          wgts_j <- xMat_train[j, ]
          if (sum(wgts_j) == 0) wgts_j <- rep(1/length(wgts_j), length(wgts_j))
          else wgts_j <- wgts_j / sum(wgts_j)
          entropic_barycenter(r_list_train, wgts_j, atoms_train, optns$lambda, optns$max_iter, optns$tol)
        })
      }
      
      # Errors
      errors <- sapply(seq_len(n_test), function(j) {
        pred_samples <- atoms_train[sample(1:K_train, size = nrow(y_test[[j]]), prob = fit_fm[[j]]), ]
        T4transport::swdist(pred_samples, y_test[[j]], p = 2, nproj = 100)$distance
      })
      outCV[f] <- mean(errors)
    }
    mean(outCV)
  }
  
  # Bandwidth ranges
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
  
  # Optimize
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

#' Conditional Wasserstein Barycenters for Multivariate Distributions
#' 
#' Implements Fréchet regression for multivariate distributional responses using Sinkhorn distance.
#' This function can handle both global and local regression methods, with automatic bandwidth
#' selection via cross-validation for the local method.
#' 
#' @param y List of length n of empirical measures; each element is an n_i x d matrix (variable sample sizes allowed)
#' @param x n x p matrix of predictors. If p = 1, a numeric vector of length n is also accepted.
#' @param xOut Optional nOut x p matrix of prediction points. If NULL, defaults to x.
#' @param optns Named list of options:
#'   \describe{
#'     \item{method}{"global" (default) or "local" — Fréchet regression type}
#'     \item{lambda}{Entropic regularization parameter (default: 0.4)}
#'     \item{grid_size}{Grid resolution per dimension (default: 101)}
#'     \item{bounds}{2 x d matrix for response domain; if NULL, uses 10%-90% quantiles per coordinate}
#'     \item{max_iter}{Maximum Dykstra iterations (default: 500)}
#'     \item{tol}{Convergence tolerance (default: 1e-6)}
#'     \item{kernel}{Kernel type for local regression: "gauss" (default), "rect", "epan", "quar", "gausvar"}
#'     \item{bw}{Bandwidth vector (length p) for local regression when method = "local"}
#'     \item{bwRange}{2 x p matrix of bandwidth ranges for selection (local method)}
#'     \item{cores}{Number of parallel cores for computation (default: max(1, detectCores()-1))}
#'     \item{verbose}{Logical; print progress (default: TRUE)}
#'   }
#' @return List containing:
#'   - pred_bary: nOut x grid_size^d matrix where each row contains probability masses for one prediction point
#'   - atoms: grid_size^d x d matrix of common discretization grid points
#'   - method: Regression method used ("global" or "local")
#'   - lambda: Entropic regularization parameter used
#'   - pred_dim: Dimension of predictors (p)
#'   - grid_size: Grid resolution per dimension
#'   - optns: Options used (including selected bandwidth for local method)
#' @references Fan, J. and Muller, H.G. (2024). Conditional Wasserstein Barycenters for Multivariate Distributions.
conditional_wb <- function(y, x, xOut = NULL, optns = list()) {
  # Input validation
  if (is.null(y) || is.null(x)) stop("Both y and x must be provided")
  if (!is.list(y)) stop("y must be a list of empirical measures")
  
  # Handle x as vector if p=1
  if (!is.matrix(x)) {
    if (is.vector(x) || is.data.frame(x)) x <- as.matrix(x)
    else stop("x must be a matrix, data.frame, or vector")
  }
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  
  n <- length(y)
  p <- ncol(x)
  if (length(y) != nrow(x)) stop("Number of observations in y and x must match")
  
  d <- ncol(y[[1]])
  
  # Handle xOut
  if (is.null(xOut)) {
    xOut <- x
    nOut <- n
  } else {
    if (!is.matrix(xOut)) {
      if (is.vector(xOut) || is.data.frame(xOut)) xOut <- as.matrix(xOut)
      else stop("xOut must be a matrix, data.frame, or vector")
    }
    if (is.null(dim(xOut))) xOut <- matrix(xOut, ncol = 1)
    nOut <- nrow(xOut)
  }
  
  # Defaults
  if (is.null(optns$method)) optns$method <- "global"
  if (is.null(optns$lambda)) optns$lambda <- 0.4
  if (is.null(optns$grid_size)) optns$grid_size <- 101
  if (is.null(optns$max_iter)) optns$max_iter <- 500
  if (is.null(optns$tol)) optns$tol <- 1e-6
  if (is.null(optns$kernel)) optns$kernel <- "gauss"
  if (is.null(optns$cores)) optns$cores <- max(1, parallel::detectCores() - 1)
  if (is.null(optns$verbose)) optns$verbose <- TRUE
  
  # Discretize
  if (optns$verbose) cat("Discretizing empirical measures...\n")
  disc <- discretize_measures(y, optns$grid_size, optns$bounds)
  atoms <- disc$atoms
  r_list <- disc$r_list
  K <- disc$K
  
  if (optns$verbose) cat("Discretized to K =", K, "atoms\n")
  
  # Weights
  if (optns$verbose) cat("Computing regression weights...\n")
  Kern <- kerFctn(optns$kernel)
  if (optns$method == "global") {
    if (optns$verbose) cat("Using global Fréchet regression\n")
    xbar <- colMeans(x)
    Sigma <- cov(x) * (n - 1) / n
    invSigma <- solve(Sigma)
    xMat <- t(apply(xOut, 1, function(xo) {
      dxo <- xo - xbar
      dxi <- sweep(x, 1, xbar, "-")
      as.vector(1 + dxi %*% (invSigma %*% dxo))
    }))
  } else {
    if (optns$verbose) cat("Using local Fréchet regression\n")
    if (is.null(optns$bw)) {
      if (optns$verbose) cat("Selecting bandwidth via CV...\n")
      optns$bw <- bwCV_fm(x, xOut, optns, y)
      if (optns$verbose) cat("Bandwidth selected:", paste(round(optns$bw, 4), collapse = " "), "\n")
    }
    if (length(optns$bw) != p) stop("bw must be vector of length p")
    
    xMat <- t(apply(xOut, 1, function(xt) {
      diffs <- sweep(x, 2, xt, "-")
      u <- sweep(diffs, 2, optns$bw, "/")
      if (p == 1) {
        Kvec <- Kern(u[, 1])
      } else {
        Kvec <- apply(u, 1, function(ui) prod(sapply(seq_len(p), function(dd) Kern(ui[dd]))))
      }
      as.vector(Kvec)
    }))
  }
  
  if (optns$verbose) cat("Weights computed\n")
  
  pred_bary <- matrix(0, nOut, K)
  
  if (optns$cores > 1) {
    # Parallel computation
    if (optns$verbose) cat(sprintf("Computing %d weighted barycenters using %d cores...\n", nOut, optns$cores))
    
    pred_bary_list <- mclapply(1:nOut, function(i) {
      wgts_i <- xMat[i, ]
      
      # Normalize
      if (sum(wgts_i) == 0) {
        wgts_i <- rep(1/n, n)
      } else {
        wgts_i <- wgts_i / sum(wgts_i)
      }
      
      # Barycenter
      entropic_barycenter(r_list, wgts_i, atoms, optns$lambda, optns$max_iter, optns$tol, FALSE)
    }, mc.cores = optns$cores)
    
    # Convert to matrix
    for (i in 1:nOut) {
      pred_bary[i, ] <- pred_bary_list[[i]]
    }
  } else {
    # Serial computation
    for (i in 1:nOut) {
      wgts_i <- xMat[i, ]
      
      # Normalize
      if (sum(wgts_i) == 0) {
        wgts_i <- rep(1/n, n)
        if (optns$verbose) warning("Zero-sum weights; using uniform")
      } else {
        wgts_i <- wgts_i / sum(wgts_i)
      }
      
      # Barycenter
      pred_bary[i, ] <- entropic_barycenter(r_list, wgts_i, atoms, optns$lambda, optns$max_iter, optns$tol, optns$verbose)
      
      if (optns$verbose && i %% 10 == 0) cat(sprintf("Computed %d/%d weighted barycenters\n", i, nOut))
    }
  }
  
  list(pred_bary = pred_bary, atoms = atoms, method = optns$method, lambda = optns$lambda,
       pred_dim = p, grid_size = optns$grid_size, optns = optns)
}
