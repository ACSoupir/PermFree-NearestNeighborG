# Neighbour counts within each radius (\<= r)

Neighbour counts within each radius (\<= r)

## Usage

``` r
nn_counts_cpp(x, y, radii)
```

## Arguments

- x, y:

  Point coordinates.

- radii:

  Radii (must be increasing for the binary-search fast path; any order
  is accepted).

## Value

An `N x length(radii)` integer matrix.
