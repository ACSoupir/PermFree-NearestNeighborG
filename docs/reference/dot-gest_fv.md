# Observed fv object for a correction

Observed fv object for a correction

## Usage

``` r
.gest_fv(X, radii = NULL, correction = "none", rmax = NULL)
```

## Arguments

- X:

  A `ppp`.

- radii:

  Radius vector (or `NULL` for the spatstat default grid).

- correction:

  One of "none", "rs", "han", "km".

- rmax:

  Optional maximum radius.

## Value

An `fv` object from \[spatstat.explore::Gest()\].
