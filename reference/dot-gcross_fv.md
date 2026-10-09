# Observed fv object for a bivariate correction

Observed fv object for a bivariate correction

## Usage

``` r
.gcross_fv(X, i, j, radii = NULL, correction = "none", rmax = NULL)
```

## Arguments

- X:

  A marked `ppp`.

- i, j:

  Mark levels.

- radii:

  Radius vector (or `NULL` for the spatstat default grid).

- correction:

  One of "none", "rs", "han", "km".

- rmax:

  Optional maximum radius.

## Value

An `fv` object from \[spatstat.explore::Gcross()\].
