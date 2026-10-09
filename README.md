# permfreeG

**Exact, permutation-free mean and variance of the nearest-neighbour G function
under complete spatial randomness.**

Instead of approximating a null distribution with thousands of random
permutations, **permfreeG** evaluates the mean and variance of the permutation
null in closed form.  A null band, pointwise z-scores and p-values are therefore
available immediately, with no Monte Carlo error.

## The result in one line

Let `X = {x_1, ..., x_N}` be fixed locations and let the marked set
S be a uniformly random n-subset of X.  Writing

* m_p(r) for the number of locations within distance r of p,
* T(B, s)(k) = C(B - k, s) / C(B, s), the probability that a uniform
  s-subset avoids k prescribed elements (the "avoidance table"),
* e_p = (n / N) (1 - T(N-1, n-1)(m_p)),
* f_pq = E[X_p X_q] as given in the vignettes,

the univariate G statistic satisfies

    E[G_S(r)]  = (1/N) sum_p [ 1 - T(N-1, n-1)(m_p(r)) ]

    Var(G_S(r)) = (1/n^2) [ sum_p e_p (1 - e_p)
                            + 2 sum_{p<q} ( f_pq - e_p e_q ) ],

and the bivariate cross statistic G_ij satisfies

    E[G_ij(r)] = (1/N) sum_p [ 1 - T(N-1, n_j)(m_p(r)) ],

with a variance of the same form using c2 = n_i (n_i - 1) / (N (N - 1)).  The
full derivations, with proofs and no code required, are in the vignettes.

Both formulas are exact for every finite sample.  The only numerical ingredient
is the avoidance table, computed once per analysis in log space.

## Installation

```{r install, eval = FALSE}
# development version
remotes::install_github("ACSoupir/PermFree-NearestNeighborG", ref = "dev",
                        build_vignettes = TRUE)
```

The package compiles a small C++ engine through Rcpp.  On macOS, make sure an
Apple (not Anaconda) compiler is first on the PATH before installing.

## Quick start

```{r example}
library(permfreeG)

# 30 points, 5 marked: exactly C(30, 5) = 142506 subsets
radii <- seq(0, 1.2, by = 0.15)
res <- exact_gest(sim_nnG, radii = radii, mark_col = "m", marks_i = "A")
as.data.frame(res)

plot(res)
```

The returned object has one row per radius with columns `r`, `theo`,
`obs`, `csr_mean`, `csr_var`, `z` and `p_value`.

For two marks, use

```{r cross, eval = FALSE}
exact_gcross(data, radii = radii, mark_col = "m", marks_i = "A",
             marks_j = "B")
```

## Edge corrections

| correction | null mean | null variance |
|---|---|---|
| `none` (raw) | exact | **exact** |
| `rs` (border / reduced sample) | ratio of exact expectations | delta method from exact moments |
| `han` (Hanisch) | mean-field using exact expected bin counts | not available; use `permute_gest()` |
| `km` (Kaplan-Meier) | mean-field using exact expected bin counts | not available; use `permute_gest()` |

Only `correction = "none"` is fully exact; the others are labelled in the
returned object and derived honestly, with their approximations stated, in the
edge-corrections vignette.

## Speed

Neighbour sets are stored as packed 64-bit bit sets, so the intersection of two
neighbourhoods is a population count over N/64 words rather than a set operation.
The complete exact variance curve for 250 points at 151 radii takes well under a
second, and `covariance = FALSE` gives an even faster independent-points
approximation.

## Vignettes

* **Derivation: univariate Gest** - full proof of the mean and variance.
* **Derivation: bivariate Gcross** - conditioning lemma, disjoint permutation,
  full proof.
* **Supplement: edge corrections** - none / rs / han / km, with derivations and a
  table of what is exact.
* **Getting started** - usage, approximations, permutation fallback and scaling.

## Data

The proprietary mIF validation data analysed in the paper is **not** part of this
repository.  A small synthetic example, `sim_nnG`, is shipped instead; it can
be regenerated with `data-raw/make-sim.R`.

The original paper scripts are preserved under `inst/legacy/` for
provenance.

## Citation

Soupir AC, Manley BJ, Peres LC, Fridley BL, Wrobel J. *Exact Expectation of
Complete Spatial Randomness for Nearest Neighbor G(r): A Scalable Alternative to
Permutations.* bioRxiv 2025.06.11.659088.
