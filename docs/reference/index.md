# Package index

## Exact permutation-free CSR

Main entry points. These compute the exact null mean and variance.

- [`exact_gest()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/exact_gest.md)
  : Permutation-free exact CSR for the univariate nearest-neighbour G
  function
- [`exact_gcross()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/exact_gcross.md)
  : Permutation-free exact CSR for bivariate cross nearest-neighbour G

## Permutation reference

Exact-by-simulation fallbacks, used to validate the closed forms.

- [`permute_gest()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/permute_gest.md)
  : Permutation null for the univariate G function
- [`permute_gcross()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/permute_gcross.md)
  : Permutation null for bivariate cross G

## Result object methods

- [`as.data.frame(`*`<nnG>`*`)`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/nnG-methods.md)
  [`print(`*`<nnG>`*`)`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/nnG-methods.md)
  [`summary(`*`<nnG>`*`)`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/nnG-methods.md)
  [`plot(`*`<nnG>`*`)`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/nnG-methods.md)
  : Methods for nnG objects

## Example data

- [`sim_nnG`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/sim_nnG.md)
  : Synthetic two-type point pattern

## Internal helpers

Not exported; used inside the package and by the test suite.

- [`avoid_table()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/avoid_table.md)
  : Avoidance table

- [`censored_neighbor_counts_cpp()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/censored_neighbor_counts_cpp.md)
  : Neighbour counts within min(radius, boundary distance)

- [`counts_within_bdist_cpp()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/counts_within_bdist_cpp.md)
  : Number of neighbours closer than the window boundary

- [`csr_cross_moments_cpp()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/csr_cross_moments_cpp.md)
  : Exact bivariate Gcross CSR mean and variance

- [`csr_moments_r()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/csr_moments_r.md)
  : Reference implementation of exact CSR moments (pure R)

- [`csr_uni_moments_cpp()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/csr_uni_moments_cpp.md)
  : Exact univariate CSR mean and variance

- [`.analytic_moments()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-analytic_moments.md)
  : Analytic CSR moments for one edge correction

- [`.check_radii()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-check_radii.md)
  : Resolve the radius grid used by spatstat

- [`.gcross_fv()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-gcross_fv.md)
  : Observed fv object for a bivariate correction

- [`.gest_fv()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-gest_fv.md)
  : Observed fv object for a correction

- [`.lapply_cores()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-lapply_cores.md)
  : Parallel lapply helper

- [`.mark_levels()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-mark_levels.md)
  : Mark levels of a point pattern

- [`.meanfield_moments()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-meanfield_moments.md)
  : Mean-field Hanisch / Kaplan-Meier CSR curve

- [`.obs_col_name()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-obs_col_name.md)
  : Column name holding the observed curve for a correction

- [`.obs_column()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-obs_column.md)
  : Pull the observed curve out of a spatstat fv object for one
  correction

- [`.prepare_ppp()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-prepare_ppp.md)
  : Coerce user data to a spatstat point pattern

- [`.subset_mark()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/dot-subset_mark.md)
  : Subset a point pattern to one mark level

- [`log_comb()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/log_comb.md)
  : Log binomial coefficient in log space

- [`new_nnG()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/new_nnG.md)
  :

  Construct an `nnG` object

- [`nn_counts_cpp()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/nn_counts_cpp.md)
  : Neighbour counts within each radius (\<= r)

- [`permfreeG`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/permfreeG-package.md)
  [`permfreeG-package`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/permfreeG-package.md)
  : permfreeG: Permutation-Free Exact Nearest-Neighbour G CSR

- [`rs_moments_cpp()`](https://acsoupir.github.io/PermFree-NearestNeighborG/reference/rs_moments_cpp.md)
  : Reduced-sample (border) CSR moments
