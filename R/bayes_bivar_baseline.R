# Bayes-Hilbert-Space Baseline for Bivariate Density Regression
# Interface matches conditional_wb from FM.R

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

#' Bayes-Hilbert-Space Baseline for Bivariate Density Regression
#'
#' Implements a Bayes-Hilbert-space-style baseline for regression from scalar
#' covariate X to bivariate density responses, using histogram discretization,
#' discrete clr transform, tensor-product spline approximation, and regression on
#' spline coefficients. Interface matches conditional_wb from FM.R.
#'
#' @param y List of length n; each element is an N x 2 numeric matrix of samples from Y_i
#' @param x n x p matrix of predictors. If p = 1, a numeric vector of length n is also accepted.
#' @param xOut Optional nOut x p matrix of prediction points. If NULL, defaults to x.
#' @param optns Named list of options:
#'   \describe{
#'     \item{method}{"global" (default) or "local" — regression type}
#'     \item{grid_size}{Grid resolution per dimension; if missing or NULL uses Sturges rule (as in Hron et al.)}
#'     \item{bounds}{2 x 2 matrix for response domain; if NULL, uses 2.5%-97.5% quantiles per coordinate}
#'     \item{eps}{Numeric; floor before log (default: 1e-12)}
#'     \item{zero_replace}{Logical; whether to apply zero replacement (default: TRUE)}
#'     \item{zero_const}{Numeric; constant for zero replacement: p0 = zero_const/(N) (default: 2/3)}
#'     \item{kx}{Integer; spline basis size in x dimension (default: 12)}
#'     \item{ky}{Integer; spline basis size in y dimension (default: 12)}
#'     \item{m_pen}{Integer vector of length 2; mgcv penalty order (default: c(1,1))}
#'     \item{lambda_spline}{Numeric or "gcv"; spline smoothing parameter (default: "gcv")}
#'     \item{lambda_grid}{Numeric vector; candidate lambdas if lambda_spline="gcv" (default: 10^seq(-4,0,len=9))}
#'     \item{kernel}{Kernel type for local regression: "gauss" (default), "rect", "epan", "quar", "gausvar"}
#'     \item{bw}{Bandwidth scalar for local regression when method = "local"}
#'     \item{verbose}{Logical; print progress (default: FALSE)}
#'   }
#'
#' @return List with components:
#'   \item{atoms}{K*L x 2 matrix of common discretization grid points (hist midpoints)}
#'   \item{predicted_densities}{List of length nOut; each element is a K x L matrix containing predicted density on grid}
#'   \item{method}{Regression method used ("global" or "local")}
#'   \item{optns}{Options used (including selected bandwidth and selected lambda if GCV)}
#'
#' @export
bayes_bivar_regression <- function(y, x, xOut = NULL, optns = list()) {
  # -------------------------
  # Input validation
  # -------------------------
  if (is.null(y) || is.null(x)) stop("Both y and x must be provided")
  if (!is.list(y)) stop("y must be a list of empirical measures")

  if (!is.matrix(x)) {
    if (is.vector(x) || is.data.frame(x)) x <- as.matrix(x)
    else stop("x must be a matrix, data.frame, or vector")
  }
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)

  n <- length(y)
  p <- ncol(x)
  if (length(y) != nrow(x)) stop("Number of observations in y and x must match")

  d <- ncol(y[[1]])
  if (d != 2) stop("This function only supports bivariate (d=2) responses")
  if (p != 1) stop("This function only supports scalar predictors (p=1)")

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

  # -------------------------
  # Defaults
  # -------------------------
  if (is.null(optns$method)) optns$method <- "global"

  # KEY: If grid_size is not provided, use Sturges rule (do NOT default to 101).
  if (!("grid_size" %in% names(optns))) optns$grid_size <- NULL

  if (is.null(optns$eps)) optns$eps <- 1e-12
  if (is.null(optns$zero_replace)) optns$zero_replace <- TRUE
  if (is.null(optns$zero_const)) optns$zero_const <- 2/3

  if (is.null(optns$kx)) optns$kx <- 12
  if (is.null(optns$ky)) optns$ky <- 12

  # Hron uses u=v=1 derivative order; mgcv analogue (in spirit) is m_pen=c(1,1).
  if (is.null(optns$m_pen)) optns$m_pen <- c(1, 1)

  # Hron selects smoothing by GCV; default to "gcv" to match their approach.
  if (is.null(optns$lambda_spline)) optns$lambda_spline <- "gcv"
  if (is.null(optns$lambda_grid)) optns$lambda_grid <- 10^seq(-4, 0, length.out = 9)

  if (is.null(optns$kernel)) optns$kernel <- "gauss"
  if (is.null(optns$verbose)) optns$verbose <- FALSE

  # Determine bounds if not provided
  if (is.null(optns$bounds)) {
    all_y <- do.call(rbind, y)
    optns$bounds <- apply(all_y, 2, quantile, probs = c(0.025, 0.975))
  }
  if (!is.matrix(optns$bounds) || nrow(optns$bounds) != 2 || ncol(optns$bounds) != 2) {
    stop("bounds must be a 2 x 2 matrix: rows = (lower, upper), cols = (dim1, dim2)")
  }

  xlim <- optns$bounds[, 1]
  ylim <- optns$bounds[, 2]

  # Default bandwidth for local regression
  if (optns$method == "local" && is.null(optns$bw)) {
    optns$bw <- 0.25 * n^(-1/5)
  }

  if (optns$verbose) cat("Building histogram grid...\n")

  # -------------------------
  # Step 1: Grid definition (hist breaks + midpoints)
  # -------------------------
  grid_obj <- make_hist_grid(xlim, ylim, optns$grid_size, y)
  gx_mid <- grid_obj$gx_mid
  gy_mid <- grid_obj$gy_mid
  x_breaks <- grid_obj$x_breaks
  y_breaks <- grid_obj$y_breaks
  K <- length(gx_mid)
  L <- length(gy_mid)
  dx <- diff(x_breaks)[1]
  dy <- diff(y_breaks)[1]

  atoms <- as.matrix(expand.grid(x = gx_mid, y = gy_mid))
  df_grid <- expand.grid(x = gx_mid, y = gy_mid)

  if (optns$verbose) cat("Building tensor-product spline basis...\n")

  # -------------------------
  # Step 2: Build basis (once)
  # -------------------------
  basis_obj <- build_basis(df_grid, optns$kx, optns$ky, optns$m_pen)
  B <- basis_obj$B
  S <- basis_obj$S
  M <- ncol(B)

  if (optns$verbose) cat("Computing discrete clr surfaces from histograms...\n")

  # -------------------------
  # Step 3: Histogram -> discrete probabilities -> discrete clr (KxL)
  # -------------------------
  Zmat <- matrix(0, nrow = n, ncol = K * L)
  for (i in seq_len(n)) {
    p_kl <- hist2_prob(y[[i]], x_breaks, y_breaks,
                       zero_replace = optns$zero_replace,
                       zero_const = optns$zero_const)
    z_kl <- clr_discrete(p_kl, eps = optns$eps)
    Zmat[i, ] <- as.vector(z_kl)  # expand.grid ordering
  }

  # -------------------------
  # Step 4: choose lambda (optional GCV)
  # -------------------------
  if (is.character(optns$lambda_spline) && tolower(optns$lambda_spline) == "gcv") {
    if (optns$verbose) cat("Selecting lambda_spline by (approx) GCV over densities...\n")
    optns$lambda_spline <- select_lambda_gcv(B, S, Zmat, optns$lambda_grid)
    if (optns$verbose) cat("Selected lambda_spline =", optns$lambda_spline, "\n")
  }

  if (optns$verbose) cat("Fitting spline coefficients...\n")

  # -------------------------
  # Step 5: Fit spline coefficients for each density
  # Uses projection to enforce mean-zero (clr) constraint in the fit:
  # minimize || P(Bc - z) ||^2 + lambda c'Sc, P = I - 11'/m.
  # -------------------------
  Cmat <- matrix(0, nrow = n, ncol = M)
  for (i in seq_len(n)) {
    Cmat[i, ] <- fit_coeff_one(B, S, Zmat[i, ], optns$lambda_spline)
  }

  if (optns$verbose) cat("Fitting regression on coefficients...\n")

  # -------------------------
  # Step 6: Regression on coefficients (global/local)
  # -------------------------
  x_train <- as.vector(x)

  if (optns$method == "global") {
    reg_params <- fit_global(x_train, Cmat)
  } else {
    Kern <- kerFctn(optns$kernel)
    reg_params <- list(bw = optns$bw, kernel = Kern, x_train = x_train, Cmat = Cmat)
  }

  if (optns$verbose) cat("Predicting densities for", nOut, "query points...\n")

  # -------------------------
  # Step 7: Predict densities for xOut
  # -------------------------
  predicted_densities <- vector("list", nOut)

  for (j in seq_len(nOut)) {
    x0 <- as.vector(xOut[j, 1])
    c_hat <- if (optns$method == "global") pred_global(x0, reg_params) else pred_local(x0, reg_params)

    # predicted clr at midpoints
    zhat_vec <- as.vector(B %*% c_hat)
    zhat <- matrix(zhat_vec, nrow = K, ncol = L)

    # enforce clr constraint on prediction (numerical guard)
    zhat <- zhat - mean(zhat)

    predicted_densities[[j]] <- inv_clr_to_density(zhat, dx, dy, eps = optns$eps)

    if (optns$verbose && j %% 10 == 0) {
      cat(sprintf("Computed %d/%d predictions\n", j, nOut))
    }
  }

  list(
    atoms = atoms,
    predicted_densities = predicted_densities,
    method = optns$method,
    optns = optns
  )
}

# ==============================================================================
# Helper Functions
# ==============================================================================

# Use Sturges if grid_size is NULL; otherwise fixed grid_size.
make_hist_grid <- function(xlim, ylim, grid_size, y_list) {
  if (is.null(grid_size)) {
    N <- nrow(y_list[[1]])
    K <- ceiling(log2(N) + 1)
    L <- K
  } else {
    K <- as.integer(grid_size)
    L <- as.integer(grid_size)
  }

  x_breaks <- seq(xlim[1], xlim[2], length.out = K + 1)
  y_breaks <- seq(ylim[1], ylim[2], length.out = L + 1)

  gx_mid <- 0.5 * (x_breaks[-1] + x_breaks[-length(x_breaks)])
  gy_mid <- 0.5 * (y_breaks[-1] + y_breaks[-length(y_breaks)])

  list(gx_mid = gx_mid, gy_mid = gy_mid, x_breaks = x_breaks, y_breaks = y_breaks)
}

# 2D histogram -> cell probability table K x L
hist2_prob <- function(samples, x_breaks, y_breaks,
                       zero_replace = TRUE, zero_const = 2/3) {
  N <- nrow(samples)

  ix <- findInterval(samples[, 1], x_breaks, rightmost.closed = TRUE)
  iy <- findInterval(samples[, 2], y_breaks, rightmost.closed = TRUE)

  K <- length(x_breaks) - 1
  L <- length(y_breaks) - 1

  ok <- (ix >= 1 & ix <= K & iy >= 1 & iy <= L)
  ix <- ix[ok]; iy <- iy[ok]

  counts <- matrix(0, nrow = K, ncol = L)
  for (t in seq_along(ix)) counts[ix[t], iy[t]] <- counts[ix[t], iy[t]] + 1

  total_count <- sum(counts)
  if (total_count == 0) {
    p <- matrix(1 / (K * L), nrow = K, ncol = L)
  } else {
    p <- counts / total_count
  }

  if (zero_replace) {
    # Practical zero replacement: replace zero cell probabilities by (2/3)/N, then renormalize.
    # (This aligns with the intended spirit of zero imputation used in their application section.)
    p0 <- (zero_const) / max(N, 1)
    p[p == 0] <- p0
    p <- p / sum(p)
  }

  p
}

# Discrete clr on KxL probability table
clr_discrete <- function(p, eps = 1e-12) {
  p <- pmax(p, eps)
  logp <- log(p)
  logp - mean(logp)
}

# Build tensor-product spline basis using mgcv
build_basis <- function(df, kx, ky, m_pen) {
  smooth_obj <- mgcv::smoothCon(
    mgcv::te(x, y, bs = c("ps", "ps"), k = c(kx, ky), m = m_pen),
    data = df,
    absorb.cons = TRUE
  )[[1]]

  B <- smooth_obj$X
  S <- smooth_obj$S[[1]] + smooth_obj$S[[2]]

  list(
    B = Matrix::Matrix(B, sparse = TRUE),
    S = Matrix::Matrix(S, sparse = TRUE),
    smooth_obj = smooth_obj
  )
}

# Fit spline coefficients with clr constraint enforced in the fit via projection P = I - 11'/m.
# We avoid forming P explicitly:
#   BtPB = BtB - (Bt u)(Bt u)' where u = 1/sqrt(m) * 1
#   BtPz = Btz - (Bt u)(u'z)
fit_coeff_one <- function(B, S, zvec, lambda) {
  m <- length(zvec)
  u <- rep(1 / sqrt(m), m)  # u'u = 1

  BtB <- Matrix::crossprod(B)
  Btz <- Matrix::crossprod(B, zvec)

  # compute Bt u (vector length M)
  Btu <- Matrix::crossprod(B, u)           # M x 1
  uz <- sum(u * zvec)                      # scalar

  # rank-1 updates
  BtPB <- BtB - (Btu %*% Matrix::t(Btu))   # M x M
  BtPz <- Btz - as.numeric(uz) * Btu       # M x 1

  lhs <- BtPB + lambda * S

  # Solve using sparse Cholesky (more stable/fast)
  fac <- Matrix::Cholesky(lhs, LDL = FALSE, perm = TRUE)
  c_hat <- Matrix::solve(fac, BtPz)

  as.vector(c_hat)
}

# Faster approximate GCV:
# For each lambda, factorize lhs once, then solve for all z's.
# We use edf = tr( (BtPB + lam S)^{-1} BtPB ) and RSS in projected space.
select_lambda_gcv <- function(B, S, Zmat, lambda_grid) {
  n <- nrow(Zmat)
  m <- ncol(Zmat)  # number of grid points

  u <- rep(1 / sqrt(m), m)

  BtB <- Matrix::crossprod(B)
  Btu <- Matrix::crossprod(B, u)
  BtPB <- BtB - (Btu %*% Matrix::t(Btu))   # constant across i

  gcv_vals <- numeric(length(lambda_grid))

  for (g in seq_along(lambda_grid)) {
    lam <- lambda_grid[g]
    lhs <- BtPB + lam * S
    fac <- Matrix::Cholesky(lhs, LDL = FALSE, perm = TRUE)

    # edf = tr( lhs^{-1} BtPB )
    A <- Matrix::solve(fac, BtPB)
    edf <- sum(Matrix::diag(A))

    # RSS in projected space
    rss <- 0
    for (i in 1:n) {
      z <- Zmat[i, ]
      # project z: Pz = z - mean(z)
      zP <- z - mean(z)

      BtPz <- Matrix::crossprod(B, zP)  # since zP already projected, no rank-1 correction needed here
      c_i <- Matrix::solve(fac, BtPz)

      zhat <- as.vector(B %*% c_i)
      zhatP <- zhat - mean(zhat)

      rss <- rss + mean((zP - zhatP)^2)
    }
    rss <- rss / n

    denom <- 1 - edf / m
    gcv_vals[g] <- if (!is.finite(denom) || denom <= 0) Inf else rss / (denom^2)
  }

  lambda_grid[which.min(gcv_vals)]
}

# Global regression on coefficients
fit_global <- function(x_train, Cmat) {
  Xdesign <- cbind(1, x_train)
  Ahat <- solve(crossprod(Xdesign), crossprod(Xdesign, Cmat))
  list(Ahat = Ahat)
}

pred_global <- function(x0, reg_params) {
  as.vector(c(1, x0) %*% reg_params$Ahat)
}

# Local linear regression without forming diag(w)
pred_local <- function(x0, reg_params) {
  bw <- reg_params$bw
  Kern <- reg_params$kernel
  x_train <- reg_params$x_train
  Cmat <- reg_params$Cmat

  u <- (x_train - x0) / bw
  w <- Kern(u)

  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) w <- rep(1, length(x_train))
  w <- w / sum(w)

  D <- cbind(1, x_train - x0)           # n x 2
  # DtWD = D' W D = (D*w)' D where W = diag(w)
  DtWD <- crossprod(D, D * w)
  # DtWC = D' W C = (D*w)' C where W = diag(w)
  DtWC <- crossprod(D, Cmat * w)

  beta <- solve(DtWD, DtWC)
  as.vector(beta[1, ])
}

# Inverse clr: exp + normalize to integrate to 1 on grid (dx*dy)
inv_clr_to_density <- function(z, dx, dy, eps = 1e-12) {
  z <- z - mean(z)
  f <- exp(z)
  f <- pmax(f, eps)
  f / (sum(f) * dx * dy)
}
