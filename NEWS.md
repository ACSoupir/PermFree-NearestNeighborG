# permfreeG 0.1.0

* First release.

* Exact closed-form mean and variance of the univariate nearest-neighbour G
  function under the mark-permutation null (`exact_gest()`).
* Exact closed-form mean and variance of the bivariate cross G function
  (`exact_gcross()`), including the i = j reduction.
* Edge corrections: raw (exact), reduced sample / border via a delta method on
  exact moments, and mean-field Hanisch and Kaplan-Meier curves.
* Packed bit-set C++ engine giving exact variance for N = 250 at 151 radii in a
  fraction of a second.
* Independent-points approximation via `covariance = FALSE`, and an
  exact-by-simulation permutation fallback (`permute_gest()`,
  `permute_gcross()`).
* Four vignettes: two complete derivations, an edge-correction supplement and a
  getting-started guide.
