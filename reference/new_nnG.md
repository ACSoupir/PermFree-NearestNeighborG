# Construct an `nnG` object

Construct an `nnG` object

## Usage

``` r
new_nnG(
  df,
  correction,
  method,
  exact,
  N,
  n_i = NA_integer_,
  n_j = NA_integer_,
  statistic = "Gest"
)
```

## Arguments

- df:

  Data frame with columns `r`, `theo`, `obs`, `csr_mean`, `csr_var`,
  `z`, `p_value`.

- correction, method:

  Settings used.

- exact:

  Logical; was the analytic result fully exact?

- N, n_i, n_j:

  Point counts.

- statistic:

  `"Gest"` or `"Gcross"`.

## Value

An object of class `nnG`.
