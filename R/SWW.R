### Sliced Wasserstein Regression for Multivariate Distributions
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

##### Helper Functions ############
# Inverse Radon transform for quantile functions (uses in-R iradon)
inverse_radon_quantiles <- function(sliced_quantiles, directions, d, grid_size = 100, xlim, ylim, cores = 1) {
  # sliced_quantiles: list of L matrices, each nOut x M (quantile functions per angle)
  # directions: L x d matrix of direction vectors
  # d: ambient data dimension (only 2D reconstruction is supported here)
  # grid_size: output image size (output_size x output_size)
  # xlim, ylim: numeric vectors of length 2 giving lower/upper bounds per axis
  
  n_out <- nrow(sliced_quantiles[[1]])
  L <- length(sliced_quantiles)
  M <- ncol(sliced_quantiles[[1]])
  x_grid <- seq(xlim[1], xlim[2], length.out = grid_size)
  y_grid <- seq(ylim[1], ylim[2], length.out = grid_size)
  
  # Compute projection angles (degrees in [0, 180)) from 2D directions
  # Only first two coordinates are used for angle definition
  theta <- {
    dirs2 <- directions[, 1:2, drop = FALSE]
    ang <- atan2(dirs2[, 2], dirs2[, 1]) * 180 / pi
    ang <- (ang %% 180)
    as.numeric(ang)
  }
  
  # UPDATED: Compute symmetric rho grid covering half-diagonal
  corners_x <- c(xlim[1], xlim[2])
  corners_y <- c(ylim[1], ylim[2])
  corner_dists <- sqrt(outer(corners_x^2, rep(1, 4)) + rep(corners_y^2, each = 4))
  half_diag <- max(corner_dists)
  
  # Pre-allocate reconstructed stack as list of matrices
  reconstructed <- vector("list", n_out)
  
  reconstructed <- mclapply(seq_len(n_out), function(i) {
    sinogram <- sapply(seq_len(L), function(j) {
      qf_row <- sliced_quantiles[[j]][i, ]
      density(qf_row, from = -half_diag, to = half_diag, bw = 'SJ', adjust = 2)$y
    })
    iradon(
      radon_image = sinogram,
      theta = theta,
      filter_name = "ramp",
      interpolation = "linear",
      x_grid = x_grid,
      y_grid = y_grid
    )
  }, mc.cores = cores)
  
  dx <- diff(x_grid)[1]
  dy <- diff(y_grid)[1]
  reconstructed <- lapply(reconstructed, function(mat) {
    mat[mat < 0] <- 0
    integral_approx <- sum(mat) * dx * dy
    if (integral_approx > 0) {
      mat / integral_approx
    } else {
      mat  # Fallback if zero
    }
  })
  
  return(list(
    densities = reconstructed,
    x_grid = x_grid,
    y_grid = y_grid
  ))
}

##### Main Function ############

swwr <- function(y, x, xOut = NULL, optns = list()) {
  ## Sliced Wasserstein Regression for Multivariate Distributions
  ## 
  ## Inputs:
  ##   y: list of length n of empirical measures; each element is an n_i x d matrix (variable sample sizes allowed)
  ##   x: n x p matrix of predictors. If p = 1, a numeric vector of length n is also accepted.
  ##   xOut: optional nOut x p matrix of prediction points. If p = 1, a numeric vector of length nOut is also accepted.
  ##   optns: named list of options
  ##     - method: "global" (default) or "local" — Fréchet regression type on slices
  ##     - directions: optional L x d matrix of slicing directions; each row must be non-zero. If provided,
  ##                   rows are normalized to unit length and L is set to nrow(directions).
  ##     - L: optional integer; number of directions to generate if 'directions' is not provided.
  ##          If both 'directions' and 'L' are missing, defaults to d * 100.
  ##     - grid_size: integer; reconstruction grid size (default 100) for 2D inverse Radon
  ##     - limits: 2 x d matrix of lower/upper bounds for reconstruction domain (only used when reconstruct = TRUE and d = 2)
  ##     - cores: integer; number of parallel workers for slicing/regression (default: max(1, detectCores()-1))
  ##     - verbose: logical; print progress (default TRUE)
  ##     - reconstruct: logical; perform 2D reconstruction from slice quantiles (default FALSE; requires d = 2)
  ##     - bw: bandwidth vector (length p) for local regression when method = "local"
  ##     - bwRange: 2 x p matrix of bandwidth ranges for selection (local method)
  ##     - kernel: kernel type for local regression: "gauss" (default), "rect", "epan", "gausvar", "quar"
  ##
  ## Outputs:
  ##   A list containing:
  ##   - sliced_quantiles: list of length L; each element is an nOut x M matrix of predicted 1D quantiles per slice
  ##   - y, x, xOut, optns: the inputs (passed through). The field optns$directions contains the L x d
  ##     unit-norm direction matrix used.
  ##   - density: optional list of length nOut with reconstructed 2D densities (when reconstruct = TRUE)
  ##   - support: optional 2 x grid_size matrix with x and y grids for reconstruction (when reconstruct = TRUE)
  
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
    # safety: ensure x is n x 1, not a length-1 scalar
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
  
  # Defaults
  if (is.null(optns$method)) optns$method <- "global"
  if (is.null(optns$reconstruct)) optns$reconstruct <- FALSE
  if (is.null(optns$grid_size)) optns$grid_size <- 100
  if (is.null(optns$cores)) optns$cores <- max(1, parallel::detectCores() - 1)
  if (is.null(optns$verbose)) optns$verbose <- TRUE
  
  if (isTRUE(optns$verbose)) {
    cat("Starting Sliced Wasserstein Regression\n")
    cat("Dimensions: n =", n, ", p =", p, ", d =", d, "\n")
    cat("Slicing directions:", optns$L, ", Spatial grid:", optns$grid_size, "\n")
  }
  
  # Step 1: Set or generate direction vectors
  if (isTRUE(optns$reconstruct) && d == 2) {
    # Special case: generate uniform directions for 2D reconstruction
    theta <- seq(0, 179, by = 1)
    directions <- cbind(cos(theta * pi / 180), sin(theta * pi / 180))
    optns$L <- nrow(directions)
  } else if (!is.null(optns$directions)) {
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
    directions <- runif_on_sphere(optns$L, d, r = 1)
  }
  
  if (isTRUE(optns$verbose)) {
    cat("Using", optns$L, "slicing directions on", d, "D unit sphere\n")
  }
  
  # Step 2: Slice data along all L directions (handles variable sample sizes internally)
  sliced_data <- slice_data(y, directions, verbose = isTRUE(optns$verbose))
  
  # Step 3: Perform Fréchet regression on each of the L slices
  if (isTRUE(optns$verbose)) {
    cat("Performing Fréchet regression on", optns$L, "slices...\n")
  }
  
  # Stack sliced data into a 3D array (M x L x n) to enable fast slice extraction
  sliced_array <- simplify2array(sliced_data)

  sliced_quantiles <- mclapply(seq_len(optns$L), function(i) {
    # Extract slice i from all observations as a matrix
    slice_y <- t(sliced_array[, i, ])  # n x M matrix
    
    # Set up options for grem/lrem
    slice_optns <- optns
    
    if (optns$method == "global") {
      result <- grem(y = slice_y, x = x, xOut = xOut, optns = slice_optns)  # from utils.R
    } else if (optns$method == "local") {
      result <- lrem(y = slice_y, x = x, xOut = xOut, optns = slice_optns)  # from utils.R
    } else {
      stop("optns$method must be either 'global' or 'local'")
    }
    
    # result is already a matrix
    return(result)
  }, mc.cores = optns$cores)
  
  if (isTRUE(optns$verbose)) {
    cat("Fréchet regression completed\n")
  }
  
  # Step 4: Reconstruct multivariate densities using inverse Radon transform (optional)
  reconstructed_densities <- NULL
  if (d == 2 && isTRUE(optns$reconstruct)) {
    if (isTRUE(optns$verbose)) {
      cat("Reconstructing multivariate densities from", optns$L, "slices...\n")
    }
    # Determine domain limits
    if (is.null(optns$limits)) {
      xr <- range(unlist(lapply(y, function(mat) mat[, 1])), na.rm = TRUE)
      yr <- range(unlist(lapply(y, function(mat) mat[, 2])), na.rm = TRUE)
      xlim <- xr
      ylim <- yr
    } else {
      if (!is.matrix(optns$limits) || nrow(optns$limits) != 2 || ncol(optns$limits) != d) {
        stop("limits must be a 2 x d matrix with lower/upper bounds per dimension")
      }
      xlim <- optns$limits[, 1]
      ylim <- optns$limits[, 2]
    }

    reconstructed_densities <- inverse_radon_quantiles(
      sliced_quantiles, directions, d, optns$grid_size, xlim = xlim, ylim = ylim, cores = optns$cores
    )
  }
  
  # Prepare output (store directions inside optns for a unified interface)
  optns$directions <- directions
  output <- list(
    sliced_quantiles = sliced_quantiles,
    y = y,
    x = x,
    xOut = xOut,
    optns = optns
  )
  if (!is.null(reconstructed_densities)) {
    output$density <- reconstructed_densities$densities
    # support as 2 x grid_size matrix: first row x-grid, second row y-grid
    output$support <- rbind(reconstructed_densities$x_grid, reconstructed_densities$y_grid)
  }
  
  if (isTRUE(optns$verbose)) {
    cat("Sliced Wasserstein Regression completed\n")
  }
  
  return(output)
}
