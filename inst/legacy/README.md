# Legacy paper scripts

These files are the original analysis and demonstration scripts that accompanied
the manuscript "Exact Expectation of Complete Spatial Randomness for Nearest
Neighbor G(r): A Scalable Alternative to Permutations"
(bioRxiv 2025.06.11.659088). They are kept verbatim for provenance and
reproducibility of the paper, and they depend on packages (for example RcppEigen)
that are not required by **permfreeG** itself.

The package reimplements and extends their functionality:

| legacy file | replacement |
|---|---|
| `functions.R` (exact_csr, exact_G_variance) | `exact_gest()`, `exact_gcross()`, and the C++ engine in `src/nn_core.cpp` |
| `compute_neighbor_counts*.cpp`, `mat_mult_eigen.cpp` | superseded by the packed-bit-set engine |
| `Demonstration Script*.R` | see the package vignettes instead |

The proprietary mIF validation data used in the paper is **not** distributed with
this package.
