# Neighbour counts within min(radius, boundary distance)

Used by the Hanisch and Kaplan-Meier moment calculations: an event is
only observed when its nearest-neighbour distance does not exceed the
point's distance to the window boundary.

## Usage

``` r
censored_neighbor_counts_cpp(x, y, bdist, radii)
```

## Arguments

- x, y:

  Point coordinates.

- bdist:

  Distance of each point to the window boundary.

- radii:

  Radii.

## Value

An `N x length(radii)` integer matrix with entry \\\\\\q \ne p : d\_{pq}
\le \min(r_k, b_p)\\\\.
