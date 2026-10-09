# Number of neighbours closer than the window boundary

For each point p, the number of other points q whose distance to p is at
most the boundary distance b_p. Used by the Kaplan-Meier / Hanisch
mean-field moment calculations.

## Usage

``` r
counts_within_bdist_cpp(x, y, bdist)
```

## Arguments

- x, y:

  Point coordinates.

- bdist:

  Distance of each point to the window boundary.

## Value

Integer vector with one entry per point.
