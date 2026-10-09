# Permutation null for the univariate G function

Exact-by-simulation reference: repeatedly draw a random subset of size
`n` from all locations and recompute the observed statistic.

## Usage

``` r
permute_gest(
  data,
  radii = NULL,
  mark_col = NULL,
  marks_i = NULL,
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

- marks_i:

  Mark level defining the subset. `NULL` uses all points.

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

\[exact_gest()\]
