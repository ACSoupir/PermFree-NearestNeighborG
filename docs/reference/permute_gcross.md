# Permutation null for bivariate cross G

Exact-by-simulation reference: permute the mark labels among all
locations, keeping the type counts fixed, and recompute \\G\_{ij}\\.

## Usage

``` r
permute_gcross(
  data,
  radii = NULL,
  mark_col = NULL,
  marks_i,
  marks_j,
  correction = c("none", "rs", "han", "km"),
  nsim = 999L,
  n_cores = 1L
)
```

## Arguments

- data:

  A `data.frame` with columns `x`, `y` (and optionally a mark column),
  or an object of class
  [`ppp`](https://rdrr.io/pkg/spatstat.geom/man/ppp.html).

- radii:

  Numeric vector of radii, starting at 0. `NULL` uses the spatstat
  default grid.

- mark_col:

  Name of the mark column when `data` is a data frame.

- marks_i, marks_j:

  Mark levels defining the two types. When they are equal the call is
  delegated to \[exact_gest()\] with `n = n_i`.

- correction:

  Edge correction: "none" (exact), "rs" (border, delta method from exact
  moments), or the mean-field "han"/"km" curves.

- nsim:

  Number of permutations.

- n_cores:

  Number of workers used for the C++ engine; currently only affects
  future parallel radius evaluation (reserved).

## Value

A list with `mean`, `var`, `r`, `nsim` and `correction`.

## See also

\[exact_gcross()\]
