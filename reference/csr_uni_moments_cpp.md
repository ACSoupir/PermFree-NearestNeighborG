# Exact univariate CSR mean and variance

Exact univariate CSR mean and variance

## Usage

``` r
csr_uni_moments_cpp(x, y, radii, n_subset, covariance = TRUE)
```

## Arguments

- x, y:

  Point coordinates of the full pattern.

- radii:

  Radii at which to evaluate the moments.

- n_subset:

  Number of marked points selected in each random subset.

- covariance:

  Include the off-diagonal covariance terms (exact variance) or only the
  marginal term (independence approximation).

## Value

A list with numeric vectors `mean` and `var`.
