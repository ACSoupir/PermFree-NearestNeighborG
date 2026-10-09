# Coerce user data to a spatstat point pattern

Coerce user data to a spatstat point pattern

## Usage

``` r
.prepare_ppp(data, mark_col = NULL)
```

## Arguments

- data:

  A `data.frame` with columns `x`, `y` (and optionally a mark column),
  or an object of class
  [`ppp`](https://rdrr.io/pkg/spatstat.geom/man/ppp.html).

- mark_col:

  Name of the column holding marks when `data` is a data frame.

## Value

A `ppp` object.
