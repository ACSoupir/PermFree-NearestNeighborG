# Analytic CSR moments for one edge correction

Analytic CSR moments for one edge correction

## Usage

``` r
.analytic_moments(
  pp,
  radii,
  n_i,
  correction = "none",
  covariance = TRUE,
  cross = FALSE,
  n_j = NULL
)
```

## Arguments

- pp:

  Full point pattern (all N locations).

- radii:

  Radii at which to evaluate the moments.

- n_i:

  Subset size (univariate) or type-i set size (bivariate).

- correction:

  One of "none", "rs", "han", "km".

- covariance:

  Include off-diagonal terms when exact.

- cross:

  Logical; bivariate Gcross?

- n_j:

  Type-j set size when `cross = TRUE`.

## Value

A list with numeric vectors `mean` and `var` (`var` is `NA` for the
mean-field han/km curves).
