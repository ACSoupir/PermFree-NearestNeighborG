# Reduced-sample (border) CSR moments

Exact first and second moments of the numerator Z and denominator Y of
the reduced-sample estimator \\\hat G\_{rs}(r) = Z/Y\\, combined into
the ratio-of-expectations mean and a first-order (delta method)
variance.

## Usage

``` r
rs_moments_cpp(x, y, bdist, radii, n_i, n_j = -1L)
```

## Arguments

- x, y:

  Point coordinates of the full pattern.

- bdist:

  Distance from each point to the window boundary.

- radii:

  Radii at which to evaluate the moments.

- n_i:

  Size of the randomly selected subset (type-i set for cross).

- n_j:

  For bivariate Gcross the size of the type-j set; use `-1` for
  univariate Gest.

## Value

A list with `mean`, `var` and the raw components `z_mean`, `z_var`,
`y_mean`, `y_var`, `cov`.
