# Reference implementation of exact CSR moments (pure R)

Computes the same quantities as \[csr_uni_moments_cpp()\] and
\[csr_cross_moments_cpp()\] using straightforward loops. Intended for
validation and exposition, not production use.

## Usage

``` r
csr_moments_r(x, y, radii, n_i = NULL, n_j = NULL, covariance = TRUE)
```

## Arguments

- x, y:

  Numeric vectors of point coordinates for the full pattern.

- radii:

  Radii at which to evaluate the moments.

- n_i:

  Subset size for univariate Gest, or type-i set size for bivariate
  Gcross.

- n_j:

  Type-j set size; `NULL` selects the univariate case.

- covariance:

  Include the off-diagonal covariance terms (exact variance) or only the
  marginal term.

## Value

A list with numeric vectors `mean` and `var`.
