# Mean-field Hanisch / Kaplan-Meier CSR curve

Uses the exact expectation of the uncensored event histogram and of the
observed-distance histogram, then assembles the corrected curve by
plug-in. The variance is not available analytically and is returned as
`NA`.

## Usage

``` r
.meanfield_moments(pp, radii, n_subset, correction)
```

## Arguments

- pp:

  Full point pattern (all N locations).

- radii:

  Radii at which to evaluate the moments.

- correction:

  One of "none", "rs", "han", "km".

## Value

A list with `mean` and an all-`NA` `var`.
