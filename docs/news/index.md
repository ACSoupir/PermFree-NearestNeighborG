# Changelog

## permfreeG 0.1.0

- First release.

- Exact closed-form mean and variance of the univariate
  nearest-neighbour G function under the mark-permutation null
  ([`exact_gest()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/exact_gest.md)).

- Exact closed-form mean and variance of the bivariate cross G function
  ([`exact_gcross()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/exact_gcross.md)),
  including the i = j reduction.

- Edge corrections: raw (exact), reduced sample / border via a delta
  method on exact moments, and mean-field Hanisch and Kaplan-Meier
  curves.

- Packed bit-set C++ engine giving exact variance for N = 250 at 151
  radii in a fraction of a second.

- Independent-points approximation via `covariance = FALSE`, and an
  exact-by-simulation permutation fallback
  ([`permute_gest()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/permute_gest.md),
  [`permute_gcross()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/permute_gcross.md)).

- Four vignettes: two complete derivations, an edge-correction
  supplement and a getting-started guide.
