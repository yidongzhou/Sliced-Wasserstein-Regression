# Inverse Radon Transform Functions (translated from Python scikit-image)

#' Convert Sinogram from Circle to Square Format
#' 
#' Helper function to convert a sinogram from circular to square format by padding
#' the original circular sinogram to fit within a square matrix. This is used
#' in the inverse Radon transform process.
#' 
#' @param sinogram A matrix representing the sinogram data (circular format)
#' @return A padded square matrix containing the original sinogram data
#' @details The function calculates the diagonal length needed for a square matrix
#' that can contain the circular sinogram, then pads the sinogram with zeros
#' to center it within the square matrix.
sinogram_circle_to_square <- function(sinogram) {
  diagonal <- ceiling(sqrt(2) * nrow(sinogram))
  old_center <- nrow(sinogram) %/% 2
  new_center <- diagonal %/% 2
  pad_before <- new_center - old_center
  # Pad the sinogram
  padded_sinogram <- matrix(0, nrow = diagonal, ncol = ncol(sinogram))
  padded_sinogram[(pad_before + 1):(pad_before + nrow(sinogram)), ] <- sinogram
  
  return(padded_sinogram)
}

#' Generate Fourier Filter for Radon Transform
#' 
#' Helper function to generate various types of Fourier filters used in the
#' inverse Radon transform process. These filters are applied in the frequency
#' domain to improve reconstruction quality.
#' 
#' @param size Integer specifying the size of the filter
#' @param filter_name Character string specifying the filter type. Options include:
#'   \itemize{
#'     \item "ramp" (default): Standard ramp filter
#'     \item "shepp-logan": Shepp-Logan filter
#'     \item "cosine": Cosine filter
#'     \item "hamming": Hamming window filter
#'     \item "hann": Hann window filter
#'     \item NULL: No filtering (returns all ones)
#'   }
#' @return A matrix of size \code{size x 1} containing the Fourier filter coefficients
#' @details The function first computes a base ramp filter using FFT, then applies
#' additional windowing functions based on the specified filter type. Different
#' filters have different frequency response characteristics that affect the
#' reconstruction quality and noise characteristics.
get_fourier_filter <- function(size, filter_name = "ramp") {
  # Create ramp filter
  n <- c(seq(1, size/2, by = 2), seq(size/2 - 1, 1, by = -2))
  f <- rep(0, size)
  f[1] <- 0.25
  f[seq(2, size, by = 2)] <- -1 / (pi * n)^2
  
  # Compute FFT
  fourier_filter <- 2 * Re(fft(f))
  
  # Apply specific filter
  if (filter_name == "ramp") {
    # Already computed
  } else if (filter_name == "shepp-logan") {
    omega <- pi * fftfreq(size)[-1]  # Remove first element to avoid division by zero
    fourier_filter[-1] <- fourier_filter[-1] * sin(omega) / omega
  } else if (filter_name == "cosine") {
    freq <- seq(0, pi, length.out = size)
    cosine_filter <- fftshift(sin(freq))
    fourier_filter <- fourier_filter * cosine_filter
  } else if (filter_name == "hamming") {
    fourier_filter <- fourier_filter * fftshift(hamming(size))
  } else if (filter_name == "hann") {
    fourier_filter <- fourier_filter * fftshift(hanning(size))
  } else if (is.null(filter_name)) {
    fourier_filter[] <- 1
  }
  
  return(matrix(fourier_filter, ncol = 1))
}

#' Generate FFT Frequency Array
#' 
#' Helper function equivalent to numpy's fftfreq. Generates the discrete
#' Fourier transform sample frequencies for a given window length.
#' 
#' @param n Integer specifying the window length
#' @param d Numeric scalar specifying the sample spacing (inverse of the sampling rate).
#'   Default is 1.0
#' @return A numeric vector of length \code{n} containing the sample frequencies
#' @details The function returns the frequencies corresponding to the discrete
#' Fourier transform bins. For even \code{n}, frequencies range from 0 to n/2-1
#' and then from -n/2 to -1. For odd \code{n}, frequencies range from 0 to (n-1)/2
#' and then from -(n-1)/2 to -1.
fftfreq <- function(n, d = 1.0) {
  if (n %% 2 == 0) {
    # Even length
    result <- c(seq(0, n/2 - 1), seq(-n/2, -1)) / (n * d)
  } else {
    # Odd length
    result <- c(seq(0, (n-1)/2), seq(-(n-1)/2, -1)) / (n * d)
  }
  return(result)
}

#' Shift Zero-Frequency Component to Center
#' 
#' Helper function equivalent to numpy's fftshift. Shifts the zero-frequency
#' component to the center of the spectrum by rearranging the elements of
#' the input array.
#' 
#' @param x A numeric vector to be shifted
#' @return A numeric vector with the zero-frequency component shifted to the center
#' @details For even-length arrays, swaps the first and second halves. For odd-length
#' arrays, swaps the first half with the second half, leaving the middle element
#' in place. This is commonly used in signal processing to center the frequency
#' spectrum around zero frequency.
fftshift <- function(x) {
  n <- length(x)
  if (n %% 2 == 0) {
    # Even length
    return(c(x[(n/2 + 1):n], x[1:(n/2)]))
  } else {
    # Odd length
    return(c(x[((n+1)/2 + 1):n], x[1:((n+1)/2)]))
  }
}

#' Generate Hamming Window
#' 
#' Helper function to generate a Hamming window, which is a type of windowing
#' function used in signal processing to reduce spectral leakage in Fourier transforms.
#' 
#' @param n Integer specifying the length of the window
#' @return A numeric vector of length \code{n} containing the Hamming window coefficients
#' @details The Hamming window is defined as w(n) = 0.54 - 0.46 * cos(2πn/(N-1))
#' for n = 0, 1, ..., N-1. It provides good frequency resolution with moderate
#' side lobe suppression. If n = 1, returns 1.
hamming <- function(n) {
  if (n == 1) return(1)
  return(0.54 - 0.46 * cos(2 * pi * seq(0, n-1) / (n-1)))
}

#' Generate Hann Window (Hanning Window)
#' 
#' Helper function to generate a Hann window (also known as Hanning window),
#' which is a type of windowing function used in signal processing to reduce
#' spectral leakage in Fourier transforms.
#' 
#' @param n Integer specifying the length of the window
#' @return A numeric vector of length \code{n} containing the Hann window coefficients
#' @details The Hann window is defined as w(n) = 0.5 * (1 - cos(2πn/(N-1)))
#' for n = 0, 1, ..., N-1. It provides good frequency resolution with better
#' side lobe suppression than the rectangular window. If n = 1, returns 1.
hanning <- function(n) {
  if (n == 1) return(1)
  return(0.5 * (1 - cos(2 * pi * seq(0, n-1) / (n-1))))
}

#' Inverse Radon Transform
#' 
#' Reconstruct an image from its Radon transform (sinogram) using filtered back-projection.
#' This function implements the inverse Radon transform with support for various filters
#' and interpolation methods.
#' 
#' @param radon_image 2D matrix representing the sinogram (Radon transform)
#' @param theta Numeric vector of projection angles in degrees. If NULL, assumes
#'   uniform angles from 0 to 180 degrees with length equal to ncol(radon_image)
#' @param filter_name Character string specifying the reconstruction filter:
#'   \itemize{
#'     \item "ramp" (default): Standard ramp filter
#'     \item "shepp-logan": Shepp-Logan filter
#'     \item "cosine": Cosine filter
#'     \item "hamming": Hamming window filter
#'     \item "hann": Hann window filter
#'     \item NULL: No filtering
#'   }
#' @param interpolation Character string specifying interpolation method:
#'   \itemize{
#'     \item "linear" (default): Linear interpolation
#'     \item "nearest": Nearest neighbor interpolation
#'     \item "cubic": Cubic spline interpolation
#'   }
#' @param circle Logical flag indicating whether to apply circular mask to output
#' @param preserve_range Logical flag indicating whether to preserve input data range
#' @param x_grid Numeric vector defining the x-coordinates of the reconstruction grid
#' @param y_grid Numeric vector defining the y-coordinates of the reconstruction grid
#' @return 2D matrix representing the reconstructed image with dimensions matching
#'   the input grids. The matrix is transposed so that rows correspond to x_grid
#'   and columns correspond to y_grid.
#' @details The function implements filtered back-projection reconstruction:
#' \enumerate{
#'   \item Applies Fourier domain filtering to the sinogram
#'   \item Performs back-projection by interpolating filtered projections
#'   \item Normalizes the result to approximate unit integral
#'   \item Optionally applies circular masking
#' }
#' The reconstruction uses physical coordinates based on the provided grids,
#' ensuring proper scaling and domain representation.
#' @references
#' \itemize{
#'   \item Kak, A.C. and Slaney, M. (2001). Principles of Computerized Tomographic Imaging.
#'   \item Bracewell, R.N. (2003). The Fourier Transform and Its Applications.
#' }
iradon <- function(radon_image, theta = NULL, 
                   filter_name = "ramp", interpolation = "linear", 
                   circle = FALSE, preserve_range = TRUE,
                   x_grid, y_grid) {  # REQUIRED: Grids determine size/domain
  # Input validation
  if (length(dim(radon_image)) != 2) {
    stop("The input image must be 2-D")
  }
  
  if (is.null(theta)) {
    theta <- seq(0, 180, length.out = ncol(radon_image))
  }
  
  angles_count <- length(theta)
  if (angles_count != ncol(radon_image)) {
    stop("The given theta does not match the number of projections in radon_image.")
  }
  
  interpolation_types <- c("linear", "nearest", "cubic")
  if (!interpolation %in% interpolation_types) {
    stop("Unknown interpolation: ", interpolation)
  }
  
  filter_types <- c("ramp", "shepp-logan", "cosine", "hamming", "hann", NULL)
  if (!filter_name %in% filter_types) {
    stop("Unknown filter: ", filter_name)
  }
  
  # Validate grids
  grid_size <- length(x_grid)
  if (length(y_grid) != grid_size) {
    stop("x_grid and y_grid must have the same length")
  }
  # Assume uniform spacing; could add diff check if needed
  dx <- diff(x_grid)[1]
  dy <- diff(y_grid)[1]
  if (is.na(dx) || dx <= 0 || is.na(dy) || dy <= 0) {
    stop("x_grid and y_grid must be strictly increasing with uniform spacing")
  }
  
  # Convert to float if needed
  if (!preserve_range) {
    radon_image <- as.double(radon_image)
  }
  storage.mode(radon_image) <- "double"
  
  img_shape <- nrow(radon_image)
  
  if (circle) {
    radon_image <- sinogram_circle_to_square(radon_image)
    img_shape <- nrow(radon_image)
  }
  
  # Resize image to next power of two (but no less than 64) for Fourier analysis
  projection_size_padded <- max(64, 2^ceiling(log2(2 * img_shape)))
  pad_width <- c(0, projection_size_padded - img_shape)
  img <- rbind(radon_image, matrix(0, nrow = pad_width[2], ncol = ncol(radon_image)))
  
  # Apply filter in Fourier domain
  fourier_filter <- get_fourier_filter(projection_size_padded, filter_name)
  
  # Apply FFT along rows (axis=0 equivalent)
  projection <- apply(img, 2, fft) * as.vector(fourier_filter)
  inv_ifft <- function(col) Re(fft(col, inverse = TRUE) / length(col))
  radon_filtered_full <- apply(projection, 2, inv_ifft)
  radon_filtered <- radon_filtered_full[1:img_shape, , drop = FALSE]
  
  # Set up physical coordinates from provided grids
  use_physical <- TRUE
  xpr_phys <- matrix(rep(x_grid, times = length(y_grid)), nrow = length(y_grid), byrow = TRUE)
  ypr_phys <- matrix(rep(y_grid, times = length(x_grid)), nrow = length(y_grid), byrow = FALSE)
  
  # Domain extents from grids
  span_x <- max(x_grid) - min(x_grid)
  span_y <- max(y_grid) - min(y_grid)
  half_diag <- sqrt((span_x / 2)^2 + (span_y / 2)^2)
  radius_phys <- min(span_x / 2, span_y / 2)
  
  detector_indices <- seq(-floor((img_shape - 1) / 2), ceiling((img_shape - 1) / 2), by = 1)
  dt <- 2 * half_diag / (img_shape - 1)
  detector_coords <- detector_indices * dt
  
  # Reconstruct image by interpolation
  reconstructed <- matrix(0, nrow = grid_size, ncol = grid_size)
  
  # Process each angle
  for (i in seq_along(theta)) {
    angle <- theta[i] * pi / 180  # Convert to radians
    col <- radon_filtered[, i]
    
    # Calculate projection coordinates
    t <- xpr_phys * cos(angle) + ypr_phys * sin(angle)
    
    # Interpolate
    if (interpolation == "linear") {
      interpolated <- approx(detector_coords, col, as.vector(t), method = "linear", 
                             yleft = 0, yright = 0, rule = 2)$y
    } else if (interpolation == "nearest") {
      interpolated <- approx(detector_coords, col, as.vector(t), method = "constant", 
                             yleft = 0, yright = 0, rule = 2)$y
    } else if (interpolation == "cubic") {
      interpolated <- spline(detector_coords, col, xout = as.vector(t), method = "natural")$y
      interpolated[is.na(interpolated)] <- 0
    }
    
    # Reshape and add to reconstruction
    reconstructed <- reconstructed + matrix(interpolated, nrow = grid_size, ncol = grid_size)
  }
  
  # Scale by normalization factor (corrected: 1 / angles_count for uniform angles over [0,180))
  reconstructed <- reconstructed * (1 / angles_count)
  
  # Apply circle mask if needed (using physical radius)
  if (circle) {
    out_reconstruction_circle <- (xpr_phys^2 + ypr_phys^2) > radius_phys^2
    reconstructed[out_reconstruction_circle] <- 0.0
  }
  
  # Normalize so that integral ≈ 1 (using actual grid spacings)
  integral_approx <- sum(reconstructed) * dx * dy
  if (integral_approx > 0) {
    reconstructed <- reconstructed / integral_approx
  }
  
  # ALWAYS: Transpose for column-major output (rows ~ x_grid, cols ~ y_grid)
  reconstructed <- t(reconstructed)
  
  return(reconstructed)
}

#' Greatest Common Divisor
#' 
#' Compute the greatest common divisor (GCD) of two integers using the Euclidean algorithm.
#' 
#' @param n Integer scalar
#' @param m Integer scalar
#' @return Integer scalar representing the greatest common divisor of \code{n} and \code{m}
#' @details The function uses the Euclidean algorithm to efficiently compute the GCD.
#' If both \code{n} and \code{m} are zero, returns 0. The function handles negative
#' inputs by taking their absolute values. Input validation ensures that both
#' arguments are integer scalars.
gcd <- function(n, m) {
  stopifnot(is.numeric(n), is.numeric(m))
  if (length(n) != 1 || floor(n) != ceiling(n) ||
      length(m) != 1 || floor(m) != ceiling(m)) {
    stop("Arguments 'n', 'm' must be integer scalars.")
  }
  if (n == 0 && m == 0) {
    return(0)
  }
  
  n <- abs(n)
  m <- abs(m)
  if (m > n) {
    t <- n
    n <- m
    m <- t
  }
  while (m > 0) {
    t <- n
    n <- m
    m <- t %% m
  }
  return(n)
}

#' Least Common Multiple
#' 
#' Compute the least common multiple (LCM) of two integers.
#' 
#' @param n Integer scalar
#' @param m Integer scalar
#' @return Integer scalar representing the least common multiple of \code{n} and \code{m}
#' @details The function computes the LCM using the relationship LCM(a,b) = |a*b|/GCD(a,b).
#' If both \code{n} and \code{m} are zero, returns 0. Input validation ensures that
#' both arguments are integer scalars.
lcm <- function(n, m) {
  stopifnot(is.numeric(n), is.numeric(m))
  if (length(n) != 1 || floor(n) != ceiling(n) ||
      length(m) != 1 || floor(m) != ceiling(m)) {
    stop("Arguments 'n', 'm' must be integer scalars.")
  }
  if (n == 0 && m == 0) {
    return(0)
  }
  
  return(n / gcd(n, m) * m)
}

#' Pairwise Least Common Multiple
#' 
#' Compute the least common multiple of a vector of integers by iteratively
#' computing the LCM of pairs of numbers.
#' 
#' @param x Numeric vector of integers
#' @return Integer scalar representing the least common multiple of all elements in \code{x}
#' @details The function computes the LCM of all elements in the input vector by
#' iteratively applying the LCM function to pairs of numbers. Zero values are
#' filtered out before computation. If the input is empty after filtering zeros,
#' returns 0. If only one non-zero element remains, returns that element.
plcm <- function(x) {
  stopifnot(is.numeric(x))
  # if (any(floor(x) != ceiling(x)) || length(x) < 2)
  #   stop("Argument 'x' must be an integer vector of length >= 2.")
  
  x <- x[x != 0]
  n <- length(x)
  if (n == 0) {
    l <- 0
  } else if (n == 1) {
    l <- x
  } else if (n == 2) {
    l <- lcm(x[1], x[2])
  } else {
    l <- lcm(x[1], x[2])
    for (i in 3:n) {
      l <- lcm(l, x[i])
    }
  }
  return(l)
}

#' Sample Uniformly on Unit Sphere
#' 
#' Generate random points uniformly distributed on the surface of a unit sphere
#' in d-dimensional space using the Box-Muller transformation.
#' 
#' @param n Integer specifying the number of points to generate
#' @param d Integer specifying the dimension of the sphere
#' @param r Numeric scalar specifying the radius of the sphere (default: 1)
#' @return A matrix of size \code{n x d} where each row represents a point on the sphere
#' @details The function generates points uniformly distributed on the surface of
#' a d-dimensional sphere by first generating multivariate normal random variables
#' and then normalizing them to unit length. This method ensures uniform distribution
#' on the sphere surface.
runif_on_sphere <- function(n, d, r = 1) {
  sims <- matrix(rnorm(n * d), nrow = n, ncol = d)
  norms <- sqrt(rowSums(sims * sims))
  r * sims / norms
}

#' Inverse Permutation
#' 
#' Compute the inverse permutation of a given permutation vector.
#' 
#' @param p Integer vector representing a permutation (indices that sort a vector)
#' @return Integer vector representing the inverse permutation
#' @details Given a permutation \code{p} such that \code{x[p]} is sorted, this function
#' returns a permutation \code{q} such that \code{q[p] = 1:length(p)}. This means
#' that if \code{y_sorted[q]} is applied to the sorted vector, it will align
#' with the original order of \code{x}.
invPerm <- function(p) {
  q <- integer(length(p))
  q[p] <- seq_along(p)
  q
}

#' Slice Multivariate Data Along Directions
#' 
#' Project multivariate data onto specified direction vectors, handling variable
#' sample sizes by extending data to a common length.
#' 
#' @param data_list List of n matrices, each n_i x d representing multivariate observations
#' @param directions L x d matrix of direction vectors for projection
#' @param verbose Logical flag to print progress messages (default: FALSE)
#' @return List of n matrices, each M x L, where M is the common length after extension
#' @details The function handles datasets with different sample sizes by extending
#' them to a common length M. If sample sizes are close (within 20% of minimum),
#' it uses the minimum size. Otherwise, it uses the least common multiple (capped
#' at 5000 for practicality). Each matrix in the output represents projections
#' of the corresponding dataset onto all L direction vectors.
slice_data <- function(data_list, directions, verbose = FALSE) {
  
  n <- length(data_list)
  L <- nrow(directions)
  d <- ncol(directions)
  
  # Handle variable sample sizes by extending to common length
  sample_sizes <- sapply(data_list, nrow)
  
  if (length(unique(sample_sizes)) == 1) {
    # All sample sizes are exactly the same - no extension needed
    M <- sample_sizes[1]
    data_extended <- data_list
  } else {
    # Sample sizes vary - need to extend
    size_range <- max(sample_sizes) - min(sample_sizes)
    size_ratio <- size_range / min(sample_sizes)
    
    if (size_ratio < 0.2) {
      # If sample sizes are close (within 10% of min), use min
      M <- min(sample_sizes)
      if (verbose) {
        cat("Sample sizes are close, using min length:", M, "\n")
      }
    } else {
      # Otherwise use least common multiple (capped for practicality)
      M <- min(plcm(sample_sizes), n * max(sample_sizes), 5000)
      if (verbose) {
        cat("Sample sizes vary significantly, using least common multiple length:", M, "\n")
      }
    }
    
    # Extend original data to common length M
    data_extended <- lapply(seq_along(data_list), function(i) {
      data_i <- data_list[[i]]
      n_i <- nrow(data_i)
      residual <- M %% n_i
      if(residual) {
        indices <- c(rep(1:n_i, each = M %/% n_i), sample(1:n_i, residual))
        data_i[indices, ]
      } else {
        rep_indices <- rep(1:n_i, each = M %/% n_i)
        data_i[rep_indices, ]
      }
    })
  }
  
  # Slice data along all directions at once -> list of n matrices M x L
  sliced_data <- lapply(seq_along(data_extended), function(j) {
    X <- data_extended[[j]]
    as.matrix(X) %*% t(directions)
  })
  
  if (verbose) {
    cat("Data slicing completed\n")
  }
  
  return(sliced_data)
}

#' @title Global Regression with Empirical Measures (REM) - Simplified
#' @description Simplified global regression for empirical measures with Euclidean predictors.
#' Assumes all empirical measures have the same length.
#' @param y A matrix of \eqn{n} empirical measures, each row is an \eqn{M}-length vector of observed values.
#' @param x An n by p matrix of predictors.
#' @param xOut An nOut by p matrix of output predictor levels. Default is \code{x}.
#' @param optns A list of options control parameters specified by \code{list(name = value)}. 
#' See `Details'.
#' @details Available control options are
#' \describe{
#' \item{lower}{The lower bound of the support of the measure. Default is \code{NULL}.}
#' \item{upper}{The upper bound of the support of the measure. Default is \code{NULL}.}
#' }
#' @return A matrix of quantile functions with dimensions nOut by M, where each row 
#' corresponds to a prediction point in xOut and each column corresponds to a quantile level.
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.G., 2023. Wasserstein Regression with Empirical Measures and Density Estimation for Sparse Data. arXiv preprint arXiv:2308.12540.}
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' }
#' @export

grem <- function(y = NULL,
                 x = NULL,
                 xOut = x,
                 optns = list()) {
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  nOut <- nrow(xOut) # number of predictions
  
  # y is a matrix (n x M)
  M <- ncol(y)
  
  # Sort each row to create quantile functions (grem expects sorted quantiles)
  y <- t(apply(y, 1, sort))  # n x M matrix with sorted rows
  
  # initialization of OSQP solver
  A <- cbind(diag(M), rep(0, M)) + cbind(rep(0, M), -diag(M))
  if (!is.null(optns$upper) &&
      !is.null(optns$lower)) {
    # if lower & upper are neither NULL
    l <- c(optns$lower, rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$upper)) {
    # if lower is NULL
    A <- A[, -1]
    l <- c(rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$lower)) {
    # if upper is NULL
    A <- A[, -ncol(A)]
    l <- c(optns$lower, rep(0, M - 1))
  } else {
    # if both lower and upper are NULL
    A <- A[, -c(1, ncol(A))]
    l <- rep(0, M - 1)
  }
  P <- diag(M)
  A <- t(A)
  q <- rep(0, M)
  u <- rep(Inf, length(l))
  model <-
    osqp::osqp(
      P = P,
      q = q,
      A = A,
      l = l,
      u = u,
      osqp::osqpSettings(max_iter = 1e05, eps_abs = 1e-05, eps_rel = 1e-05, verbose = FALSE)
    )
  
  xMean <- colMeans(x)
  xc <- sweep(x, 2, xMean, FUN = "-")
  invVa <- solve(var(x) * (n - 1) / n)
  wc <- xc %*% invVa  # n x p
  
  # Only compute predictions for xOut
  qp <- matrix(nrow = nOut, ncol = M)
  for (i in 1:nOut) {
    delta <- xOut[i, ] - xMean            # 1 x p
    w <- as.vector(1 + wc %*% delta)      # n
    # column-weighted means: (t(y) %*% w) / sum(w)
    qNew <- as.vector((t(y) %*% w) / sum(w))
    if (any(w < 0)) {
      # if negative weights exist
      model$Update(q = -qNew)
      qNew <- sort(model$Solve()$x)
    }
    if (!is.null(optns$upper)) {
      qNew <- pmin(qNew, optns$upper)
    }
    if (!is.null(optns$lower)) {
      qNew <- pmax(qNew, optns$lower)
    }
    qp[i, ] <- qNew
  }
  
  # Return qp matrix directly
  return(qp)
}

#' Kernel Function Factory
#' 
#' Create kernel functions for local regression based on the specified kernel type.
#' 
#' @param kernel_type Character string specifying the kernel type. Options include:
#'   \itemize{
#'     \item "gauss": Gaussian kernel (default)
#'     \item "rect": Rectangular kernel
#'     \item "epan": Epanechnikov kernel
#'     \item "gausvar": Gaussian variance kernel
#'     \item "quar": Quartic kernel
#'   }
#' @return A function that takes a numeric vector and returns kernel weights
#' @details Different kernel types have different shapes and properties:
#' \itemize{
#'   \item Gaussian: Smooth, infinite support
#'   \item Rectangular: Constant within support, zero outside
#'   \item Epanechnikov: Parabolic shape, optimal for density estimation
#'   \item Gaussian variance: Modified Gaussian with variance adjustment
#'   \item Quartic: Fourth-order polynomial kernel
#' }
kerFctn <- function(kernel_type){
  if (kernel_type=='gauss'){
    ker <- function(x){
      dnorm(x) #exp(-x^2 / 2) / sqrt(2*pi)
    }
  } else if(kernel_type=='rect'){
    ker <- function(x){
      as.numeric((x<=1) & (x>=-1))
    }
  } else if(kernel_type=='epan'){
    ker <- function(x){
      n <- 1
      (2*n+1) / (4*n) * (1-x^(2*n)) * (abs(x)<=1)
    }
  } else if(kernel_type=='gausvar'){
    ker <- function(x) {
      dnorm(x)*(1.25-0.25*x^2)
    }
  } else if(kernel_type=='quar'){
    ker <- function(x) {
      (15/16)*(1-x^2)^2 * (abs(x)<=1)
    }
  } else {
    stop('Unavailable kernel')
  }
  return(ker)
}

#' Set Bandwidth Range for Kernel Regression
#' 
#' Determine appropriate bandwidth range for kernel regression based on input data
#' characteristics and kernel type.
#' 
#' @param xin Numeric vector of input predictor values
#' @param xout Numeric vector of output predictor values (prediction points)
#' @param kernel_type Character string specifying the kernel type
#' @return List with components:
#'   \item{min}{Minimum bandwidth value}
#'   \item{max}{Maximum bandwidth value}
#' @details The function calculates bandwidth bounds based on data spacing and
#' kernel-specific scaling factors. For Gaussian kernels, a scaling factor of 3
#' is applied, while Gaussian variance kernels use a factor of 2.5. The minimum
#' bandwidth ensures adequate smoothing, while the maximum prevents over-smoothing.
SetBwRange <- function(xin, xout, kernel_type) {
  xinSt <- unique(sort(xin))
  bw.min <- max(diff(xinSt), xinSt[2] - min(xout), max(xout) -
                  xinSt[length(xinSt) - 1]) * 1.1 / (ifelse(kernel_type == "gauss", 3, 1) *
                                                       ifelse(kernel_type == "gausvar", 2.5, 1))
  bw.max <- diff(range(xin)) / 3
  if (bw.max < bw.min) {
    if (bw.min > bw.max * 3 / 2) {
      # warning("Data is too sparse.")
      bw.max <- bw.min * 1.01
    } else {
      bw.max <- bw.max * 3 / 2
    }
  }
  return(list(min = bw.min, max = bw.max))
}

#' Bandwidth Selection via Cross-Validation for Local REM
#' 
#' Select optimal bandwidth for local regression with empirical measures using
#' k-fold cross-validation.
#' 
#' @param xin n x p matrix of input predictors
#' @param qin n x M matrix of quantile functions (empirical measures)
#' @param xout nOut x p matrix of output predictor levels (default: xin)
#' @param optns List of options including kernel type and bandwidth constraints
#' @return Numeric vector of length p containing optimal bandwidth values
#' @details The function performs k-fold cross-validation (k=10 for n>30, leave-one-out
#' otherwise) to minimize the mean squared error between predicted and observed
#' quantile functions. For multivariate predictors (p>1), bandwidth is optimized
#' jointly across all dimensions. The optimization uses L-BFGS-B method for
#' constrained optimization within the specified bandwidth range.
bwCV_lrem <- function(xin, qin, xout, optns) {
  if(is.null(xout)) {
    xout <- xin
  }
  n <- nrow(xin)
  p <- ncol(xin)
  # initialization of OSQP solver
  M <- ncol(qin)
  A <- cbind(diag(M), rep(0, M)) + cbind(rep(0, M), -diag(M))
  if (!is.null(optns$upper) &&
      !is.null(optns$lower)) {
    # if lower & upper are neither NULL
    l <- c(optns$lower, rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$upper)) {
    # if lower is NULL
    A <- A[, -1]
    l <- c(rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$lower)) {
    # if upper is NULL
    A <- A[, -ncol(A)]
    l <- c(optns$lower, rep(0, M - 1))
  } else {
    # if both lower and upper are NULL
    A <- A[, -c(1, ncol(A))]
    l <- rep(0, M - 1)
  }
  # P <- as(diag(M), "sparseMatrix")
  # A <- as(t(A), "sparseMatrix")
  P <- diag(M)
  A <- t(A)
  q <- rep(0, M)
  u <- rep(Inf, length(l))
  model <-
    osqp::osqp(
      P = P,
      q = q,
      A = A,
      l = l,
      u = u,
      osqp::osqpSettings(verbose = FALSE)
    )
  
  # select kernel
  Kern <- kerFctn(optns$kernel)
  K <- function(x, h) {
    k <- 1
    for (i in 1:p) {
      k <- k * Kern(x[i] / h[i])
    }
    return(as.numeric(k))
  }
  
  # k-fold
  objFctn <- function(h) {
    numFolds <- ifelse(n > 30, 10, n)# leave-one-out or 10-fold cross-validation
    folds <- sample(c(rep.int(1:numFolds, n%/%numFolds), seq_len(n%%numFolds)))
    
    cv <- 0
    for (foldidx in seq_len(numFolds)) {
      # nn by M
      testidx <- which(folds == foldidx)
      for (j in testidx) {
        a <- xin[j, ]
        if (p > 1) {
          mu1 <-
            rowMeans(apply(xin[-testidx, ], 1, function(xi)
              K(xi - a, h) * (xi - a)))
          mu2 <-
            matrix(rowMeans(apply(xin[-testidx, ], 1, function(xi)
              K(xi - a, h) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <-
            mean(sapply(xin[-testidx, ], function(xi)
              K(xi - a, h) * (xi - a)))
          mu2 <-
            mean(sapply(xin[-testidx, ], function(xi)
              K(xi - a, h) * (xi - a)^2))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(as.matrix(xin[-testidx, ]), 1, function(xi) {
          K(xi - a, h) * (1 - wc %*% (xi - a))
        }) # weight
        qNew <- apply(qin[-testidx,], 2, weighted.mean, w) # N
        if (any(w < 0)) {
          # if negative weights exist
          model$Update(q = -qNew)
          cv <- cv + sum((qin[j, ] - sort(model$Solve()$x))^2) / (n * M)
        } else {
          cv <- cv + sum((qin[j, ] - qNew)^2) / (n * M)
        }
      }
    }
    cv
  }
  
  if (p == 1) {
    aux <-
      SetBwRange(xin = xin[, 1],
                 xout = xout[, 1],
                 kernel_type = optns$kernel)
    bwRange <- matrix(c(aux$min, aux$max), nrow = 2, ncol = 1)
  } else {
    aux <-
      SetBwRange(xin = xin[, 1],
                 xout = xout[, 1],
                 kernel_type = optns$kernel)
    aux2 <-
      SetBwRange(xin = xin[, 2],
                 xout = xout[, 2],
                 kernel_type = optns$kernel)
    bwRange <-
      as.matrix(cbind(c(aux$min, aux$max), c(aux2$min, aux2$max)))
  }
  if (!is.null(optns$bw)) {
    if (p == 1) {
      if (min(optns$bw) < min(bwRange)) {
        message("Minimum bandwidth is too small and has been reset.")
      } else {
        bwRange[1, 1] <- min(optns$bw)
      }
      if (max(optns$bw) > min(bwRange)) {
        bwRange[2, 1] <- max(optns$bw)
      } else {
        message("Maximum bandwidth is too small and has been reset.")
      }
    } else {
      # Check for first dimension of the predictor
      if (min(optns$bw[, 1]) < min(bwRange[, 1])) {
        message("Minimum bandwidth of first predictor dimension is too small and has been reset.")
      } else {
        bwRange[1, 1] <- min(optns$bw[, 1])
      }
      if (max(optns$bw[, 1]) > min(bwRange[, 1])) {
        bwRange[2, 1] <- max(optns$bw[, 1])
      } else {
        message("Maximum bandwidth of first predictor dimension is too small and has been reset.")
      }
      # Check for second dimension of the predictor
      if (min(optns$bw[, 2]) < min(bwRange[, 2])) {
        message("Minimum bandwidth of second predictor dimension is too small and has been reset.")
      } else {
        bwRange[1, 2] <- min(optns$bw[, 2])
      }
      if (max(optns$bw[, 2]) > min(bwRange[, 2])) {
        bwRange[2, 2] <- max(optns$bw[, 2])
      } else {
        message("Maximum bandwidth of second predictor dimension is too small and has been reset.")
      }
    }
  }
  if (p == 1) {
    res <- optimize(f = objFctn, interval = bwRange[, 1])$minimum
  } else {
    res <-
      optim(
        par = colMeans(bwRange),
        fn = objFctn,
        lower = bwRange[1, ],
        upper = bwRange[2, ],
        method = "L-BFGS-B"
      )$par
  }
  res
}

#' @title Local Regression with Empirical Measures (REM) - Simplified
#' @description Simplified local regression for empirical measures with Euclidean predictors.
#' Assumes all empirical measures have the same length.
#' @param y A matrix of \eqn{n} empirical measures, each row is an \eqn{M}-length vector of observed values.
#' @param x An n by p matrix of predictors.
#' @param xOut An nOut by p matrix of output predictor levels. Default is \code{x}.
#' @param optns A list of options control parameters specified by \code{list(name = value)}. 
#' See `Details'.
#' @details Available control options are
#' \describe{
#' \item{bw}{A vector of length p used as the bandwidth for local REM. If not 
#' specified, the bandwidth will be selected using cross-validation.}
#' \item{bwRange}{A 2 by p matrix whose columns contain the bandwidth selection 
#' range for each corresponding dimension of the predictor \code{x} for the case 
#' when \code{bw} is \code{NULL}. Default is \code{NULL} and is automatically 
#' chosen by a data-adaptive method.}
#' \item{kernel}{A character holding the type of kernel functions for local REM; 
#' \code{"rect"}, \code{"gauss"}, \code{"epan"}, \code{"gausvar"}, \code{"quar"} 
#' - default: \code{"gauss"}.}
#' \item{lower}{The lower bound of the support of the measure. Default is \code{NULL}.}
#' \item{upper}{The upper bound of the support of the measure. Default is \code{NULL}.}
#' }
#' @return A matrix of quantile functions with dimensions nOut by M, where each row 
#' corresponds to a prediction point in xOut and each column corresponds to a quantile level.
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.G., 2023. Wasserstein Regression with Empirical Measures and Density Estimation for Sparse Data. arXiv preprint arXiv:2308.12540.}
#' \item \cite{Petersen, A. and Müller, H.-G. (2019). Fréchet regression for random objects with Euclidean predictors. The Annals of Statistics, 47(2), 691--719.}
#' }
#' @export

lrem <- function(y = NULL,
                 x = NULL,
                 xOut = x,
                 optns = list()) {
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  nOut <- nrow(xOut) # number of predictions
  
  if (!is.null(optns$bw)) {
    if (sum(optns$bw <= 0) > 0) {
      stop("bandwidth must be positive")
    }
    if (length(optns$bw) != p) {
      stop("dimension of bandwidth does not agree with x")
    }
  }
  if (!is.null(optns$bwRange)) {
    if (!is.matrix(optns$bwRange) && !is.vector(optns$bwRange)) {
      stop("bwRange must be a matrix or vector")
    }
    if (is.vector(optns$bwRange)) {
      optns$bwRange <- matrix(optns$bwRange, length(optns$bwRange))
      if (ncol(x) > 1) {
        stop("bwRange must be a matrix")
      } else {
        if (nrow(optns$bwRange) != 2) {
          stop("bwRange must have the lower and upper bound for the bandwidth range")
        }
      }
    } else {
      if (ncol(optns$bwRange) != ncol(x)) {
        stop("bwRange must have the same number of columns as x")
      }
      if (nrow(optns$bwRange) != 2) {
        stop("bwRange must have two rows")
      }
    }
  }
  
  if (is.null(optns$kernel)) {
    optns$kernel <- "gauss"
  }
  
  # y is a matrix (n x M)
  M <- ncol(y)
  
  # Sort each row to create quantile functions (lrem expects sorted quantiles)
  y <- t(apply(y, 1, sort))  # n x M matrix with sorted rows
  
  # initialization of OSQP solver
  A <- cbind(diag(M), rep(0, M)) + cbind(rep(0, M), -diag(M))
  if (!is.null(optns$upper) &&
      !is.null(optns$lower)) {
    # if lower & upper are neither NULL
    l <- c(optns$lower, rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$upper)) {
    # if lower is NULL
    A <- A[, -1]
    l <- c(rep(0, M - 1), -optns$upper)
  } else if (!is.null(optns$lower)) {
    # if upper is NULL
    A <- A[, -ncol(A)]
    l <- c(optns$lower, rep(0, M - 1))
  } else {
    # if both lower and upper are NULL
    A <- A[, -c(1, ncol(A))]
    l <- rep(0, M - 1)
  }
  P <- diag(M)
  A <- t(A)
  q <- rep(0, M)
  u <- rep(Inf, length(l))
  model <-
    osqp::osqp(
      P = P,
      q = q,
      A = A,
      l = l,
      u = u,
      osqp::osqpSettings(max_iter = 1e05, eps_abs = 1e-05, eps_rel = 1e-05, verbose = FALSE)
    )
  
  # select kernel
  Kern <- kerFctn(optns$kernel)
  
  if (is.null(optns$bw)) {
    optns$bw <- bwCV_lrem(
      xin = x,
      qin = y,
      xout = xOut,
      optns = optns
    )
  } else {
    if (ncol(x) == 1) {
      if (optns$bw[1] < max(diff(sort(x[, 1]))) &&
          !is.null(optns$kernel)) {
        if (optns$kernel %in% c("rect", "quar", "epan")) {
          warning("optns$bw was set too small and is reset to be chosen by CV.")
          optns$bw <-
            bwCV_lrem(
              xin = x,
              qin = y,
              xout = xOut,
              optns = optns
            )
        }
      }
    } else {
      if (optns$bw[1] < max(diff(sort(x[, 1]))) &&
          optns$bw[2] < max(diff(sort(x[, 2]))) && !is.null(optns$kernel)) {
        if (optns$kernel %in% c("rect", "quar", "epan")) {
          warning("optns$bw was set too small and is reset to be chosen by CV.")
          optns$bw <-
            bwCV_lrem(
              xin = x,
              qin = y,
              xout = xOut,
              optns = optns
            )
        }
      }
    }
  }
  
  # Only compute predictions for xOut
  qp <- matrix(nrow = nOut, ncol = M)
  for (i in 1:nOut) {
    a <- xOut[i, ]
    diffs <- sweep(x, 2, a, "-")
    u <- sweep(diffs, 2, optns$bw, "/")
    ker_mat <- sapply(1:p, function(d) Kern(u[, d]))
    Kvec <- apply(ker_mat, 1, prod)
    k_diffs <- diffs * Kvec
    mu1 <- colMeans(k_diffs)
    mu2_raw <- t(diffs) %*% k_diffs
    mu2 <- mu2_raw / n
    beta <- solve(mu2, mu1)
    linear_terms <- as.vector(diffs %*% beta)
    w <- Kvec * (1 - linear_terms)
    qNew <- as.vector((t(y) %*% w) / sum(w))   # M
    if (any(w < 0)) {
      # if negative weights exist
      model$Update(q = -qNew)
      qNew <- sort(model$Solve()$x)
    }
    if (!is.null(optns$upper)) {
      qNew <- pmin(qNew, optns$upper)
    }
    if (!is.null(optns$lower)) {
      qNew <- pmax(qNew, optns$lower)
    }
    qp[i, ] <- qNew
  }
  
  # Return qp matrix directly
  return(qp)
}
