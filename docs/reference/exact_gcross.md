# Permutation-free exact CSR for bivariate cross nearest-neighbour G

Computes the observed \\G\_{ij}(r)\\ together with its exact mean and
variance under the disjoint label-permutation null: a random type-i set
of size \\n_i\\ and, from the remaining locations, a random type-j set
of size \\n_j\\.

## Usage

``` r
exact_gcross(
  data,
  radii = NULL,
  mark_col = NULL,
  marks_i,
  marks_j,
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

- marks_i, marks_j:

  Mark levels defining the two types. When they are equal the call is
  delegated to \[exact_gest()\] with `n = n_i`.

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

An object of class `nnG` (see \[exact_gest()\]).

## See also

\[exact_gest()\], \[permute_gcross()\]

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(60), y = runif(60),
                m = sample(c("A", "B", "C"), 60, replace = TRUE))
res <- exact_gcross(d, radii = seq(0, 0.5, by = 0.1), mark_col = "m",
                    marks_i = "A", marks_j = "B")
as.data.frame(res)
#>     r      theo       obs  csr_mean      csr_var          z    p_value
#> 1 0.0 0.0000000 0.0000000 0.0000000 0.000000e+00         NA         NA
#> 2 0.1 0.5265373 0.3333333 0.4244506 1.274591e-02 -0.8070773 0.41962195
#> 3 0.2 0.9497492 0.6666667 0.8590137 9.088537e-03 -2.0176159 0.04363128
#> 4 0.3 0.9988044 1.0000000 0.9805102 1.615970e-03  0.4848305 0.62779656
#> 5 0.4 0.9999936 1.0000000 0.9981265 1.541219e-04  0.1509149 0.88004287
#> 6 0.5 1.0000000 1.0000000 0.9998224 1.639592e-05  0.0438521 0.96502229
```
