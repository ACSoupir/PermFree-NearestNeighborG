# Permutation-free exact CSR for the univariate nearest-neighbour G function

Computes the observed \\G(r)\\ for a marked subset together with the
exact mean and variance of that same statistic under random relabelling
(a uniform subset of the observed locations), removing any need to
enumerate permutations.

## Usage

``` r
exact_gest(
  data,
  radii = NULL,
  mark_col = NULL,
  marks_i = NULL,
  correction = c("none", "rs", "han", "km"),
  method = c("exact", "approx"),
  covariance = TRUE,
  n_cores = 1L,
  rmax = NULL
)
```

## Arguments

- data:

  A `data.frame` with columns `x`, `y` (and optionally a mark column),
  or an object of class
  [`ppp`](https://rdrr.io/pkg/spatstat.geom/man/ppp.html).

- radii:

  Numeric vector of radii, starting at 0. `NULL` uses the spatstat
  default grid.

- mark_col:

  Name of the mark column when `data` is a data frame.

- marks_i:

  Mark level defining the subset. `NULL` uses all points.

- correction:

  Edge correction: "none" (exact), "rs" (border, delta method from exact
  moments), or the mean-field "han"/"km" curves.

- method:

  `"exact"` (include the covariance terms) or `"approx"` (marginal
  variance only).

- covariance:

  Include off-diagonal covariance terms. Ignored when
  `method = "approx"`.

- n_cores:

  Number of workers used for the C++ engine; currently only affects
  future parallel radius evaluation (reserved).

- rmax:

  Optional maximum radius passed to spatstat.

## Value

An object of class `nnG`: a data frame with columns `r`, `theo`, `obs`,
`csr_mean`, `csr_var`, `z` and `p_value`.

## See also

\[exact_gcross()\], \[permute_gest()\]

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(40), y = runif(40),
                m = rep(c("A", "B"), each = 20))
res <- exact_gest(d, radii = seq(0, 0.5, by = 0.1), mark_col = "m",
                  marks_i = "A")
as.data.frame(res)
#>     r      theo  obs  csr_mean      csr_var          z    p_value
#> 1 0.0 0.0000000 0.00 0.0000000 0.000000e+00         NA         NA
#> 2 0.1 0.5942823 0.65 0.4214137 1.517436e-02 1.85564556 0.06350409
#> 3 0.2 0.9729046 0.90 0.8929692 4.275720e-03 0.10752279 0.91437424
#> 4 0.3 0.9997021 1.00 0.9835419 6.860984e-04 0.62832871 0.52978862
#> 5 0.4 0.9999995 1.00 0.9966441 1.569003e-04 0.26791683 0.78876334
#> 6 0.5 1.0000000 1.00 0.9999411 2.940538e-06 0.03433324 0.97261142
```
