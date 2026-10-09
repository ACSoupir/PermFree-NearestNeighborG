# Exact bivariate Gcross CSR mean and variance

Disjoint label-permutation null: type-i set of size `n_i`, type-j set of
size `n_j` drawn without replacement from the remaining locations.

## Usage

``` r
csr_cross_moments_cpp(x, y, radii, n_i, n_j, covariance = TRUE)
```

## Arguments

- x, y:

  Point coordinates of the full pattern.

- radii:

  Radii at which to evaluate the moments.

- n_i:

  Number of type-i points in each permutation.

- n_j:

  Number of type-j points in each permutation.

- covariance:

  Include the off-diagonal covariance terms (exact variance) or only the
  marginal term.

## Value

A list with numeric vectors `mean` and `var`.
