# Sliced Wasserstein Regression

Code for the JMLR paper **"Sliced Wasserstein Regression"** by Han Chen,
Yidong Zhou, and Hans-Georg Muller.

This repository implements regression methods for **multivariate distributional
responses** with Euclidean predictors. Each response is an empirical sample from
a multivariate probability distribution, and the goal is to estimate how the
entire conditional distribution changes with predictors.

The paper develops two slicing-based Frechet regression approaches:

- **SWW regression** (`swwr`): runs univariate Wasserstein/Frechet regression
  independently on each projection direction, then optionally reconstructs a
  bivariate density by a regularized inverse Radon transform. This is the faster
  and typically more accurate method in the paper's simulations.
- **SAW regression** (`sawr`): directly fits a sliced Wasserstein barycenter in
  the original multivariate response space using a gradient-based optimization.

Both methods support global and local Frechet regression weights.

## Repository Layout

```text
R/
  load_swr.R                 # convenience loader for all project functions
  SWW.R                      # slice-wise Wasserstein regression: swwr()
  SAW.R                      # slice-averaged/direct regression: sawr()
  FM.R                       # Fan-Muller conditional W2 barycenter baseline
  bayes_bivar_baseline.R     # Bayes-Hilbert bivariate density baseline
  utils.R                    # shared utilities, quantile regression, inverse Radon

scripts/
  sim_global.R               # global SWW/SAW simulations
  sim_local.R                # local SWW/SAW simulations
  sim_slicing_directions.R   # sensitivity to the number of slicing directions
  sim_fm.R                   # Fan-Muller baseline simulations
  sim_bayes.R                # Bayes-Hilbert baseline simulations
  runtime.R                  # runtime comparison
  visual_comparison_saw_sww.R

data/
  *.rds                      # saved simulation outputs
```

## Installation

This is a source-code repository rather than an R package. Clone it, then load
the code from the repository root:

```r
source("R/load_swr.R")
```

The current scripts use these R packages:

```r
install.packages(c(
  "MASS", "T4transport", "caret", "osqp", "ks", "ggplot2",
  "ggdensity", "mvtnorm", "scico", "mgcv", "Matrix"
))
```

`parallel` is part of base R.

## Data Format

The main functions expect:

- `y`: a list of length `n`; each element is an `N_i x d` matrix containing an
  empirical sample from the `i`th multivariate response distribution.
- `x`: an `n x p` predictor matrix. A numeric vector is accepted for scalar
  predictors.
- `xOut`: optional prediction points. If omitted, predictions are returned at
  the observed predictors.
- `optns`: a named list of method options.

## Quick Start

```r
source("R/load_swr.R")
library(MASS)

set.seed(1)
n <- 30
d <- 2
N <- 100

x <- runif(n, -0.5, 0.5)
y <- lapply(x, function(xi) {
  MASS::mvrnorm(
    n = N,
    mu = c(xi, xi),
    Sigma = diag(rep(xi + 1, d))
  )
})

xOut <- seq(-0.5, 0.5, length.out = 25)

fit_sww <- swwr(
  y = y,
  x = x,
  xOut = xOut,
  optns = list(method = "global", L = 200, cores = 1, verbose = FALSE)
)

length(fit_sww$sliced_quantiles)
fit_sww$optns$directions[1:5, ]
```

For bivariate responses, `swwr` can also reconstruct densities:

```r
fit_sww_density <- swwr(
  y = y,
  x = x,
  xOut = c(0),
  optns = list(
    method = "global",
    reconstruct = TRUE,
    grid_size = 80,
    limits = matrix(c(-4, 4, -4, 4), nrow = 2),
    cores = 1,
    verbose = FALSE
  )
)

image(fit_sww_density$support[1, ], fit_sww_density$support[2, ],
      fit_sww_density$density[[1]])
```

The direct SAW fit has a similar interface:

```r
fit_saw <- sawr(
  y = y,
  x = x,
  xOut = xOut,
  optns = list(method = "global", L = 200, cores = 1, verbose = FALSE)
)

length(fit_saw$predicted_samples)
```

For local regression, use `method = "local"` and provide a bandwidth, or let the
function select one by cross-validation:

```r
fit_local <- swwr(
  y = y,
  x = x,
  xOut = xOut,
  optns = list(method = "local", bw = n^(-0.2) / 4, cores = 1)
)
```

## Reproducing Simulations

The simulation scripts are designed to match the paper's Monte Carlo studies.
They are computationally intensive by default: most use `num_runs <- 100`, sample
sizes `n = 50, 100, 200`, dimensions `d = 2, 5`, and up to 60 worker cores. For a
quick local check, reduce `num_runs`, `dims`, `n_vals`, and `max_cores` at the top
of the script.

Run from the repository root:

```sh
Rscript scripts/sim_global.R
Rscript scripts/sim_local.R
Rscript scripts/sim_slicing_directions.R
Rscript scripts/runtime.R
```

Or from inside `scripts/`:

```sh
Rscript sim_global.R
```

All scripts now discover the project root automatically. Simulation outputs are
written to `data/`. The visual comparison script writes PDFs to `figures/`.

## Outputs

`swwr()` returns a list with:

- `sliced_quantiles`: one `nOut x M` fitted quantile matrix per slicing direction.
- `optns$directions`: the unit directions used for slicing.
- `density` and `support`: optional bivariate reconstructions when
  `reconstruct = TRUE`.

`sawr()` returns:

- `predicted_samples`: one fitted empirical sample per prediction point.
- `predicted_densities`: optional kernel density estimates when `density = TRUE`.
- `optns$directions`: the directions used in the sliced Wasserstein objective.

## Method Summary

The Radon transform maps a multivariate distribution into a collection of
one-dimensional projected distributions. This makes Wasserstein regression more
tractable because one-dimensional Wasserstein geometry can be represented through
quantile functions.

SWW uses this idea slice by slice: for each direction, project all empirical
measures, fit a univariate distributional regression, and collect the fitted
quantiles. In two dimensions, the fitted slices can be reconstructed into a
density by a regularized inverse Radon transform.

SAW keeps the response in the multivariate distribution space and minimizes a
slice-averaged Wasserstein objective directly. It is closer to a sliced
Wasserstein barycenter computation and is therefore more expensive.

The paper studies global and local versions of both estimators, proves
asymptotic guarantees, and evaluates them against Fan-Muller conditional
Wasserstein barycenters and a Bayes-Hilbert bivariate density baseline.

## Citation

```bibtex
@article{chen2026sliced,
  title = {Sliced Wasserstein Regression},
  author = {Chen, Han and Zhou, Yidong and Muller, Hans-Georg},
  journal = {Journal of Machine Learning Research},
  volume = {27},
  pages = {1--67},
  year = {2026}
}
```

## License

This code is released under the MIT License.
