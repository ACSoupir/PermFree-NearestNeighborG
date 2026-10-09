# Pull the observed curve out of a spatstat fv object for one correction

Pull the observed curve out of a spatstat fv object for one correction

## Usage

``` r
.obs_column(obs, correction)
```

## Arguments

- obs:

  An `fv` object returned by Gest/Gcross.

- correction:

  One of "none", "rs", "han", "km".

## Value

Numeric vector of observed values.
