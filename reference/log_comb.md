# Log binomial coefficient in log space

`log_comb(n, k) = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)`.
Used to avoid overflow for large `n`.

## Usage

``` r
log_comb(n, k)
```

## Arguments

- n, k:

  Numeric vectors (recycled).

## Value

Numeric vector of log-binomial coefficients.
