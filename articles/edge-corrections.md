# Supplement: edge corrections for the permutation-free G function

This vignette supplements the two derivation vignettes. It defines each
edge correction supported by the package, states exactly which part of
its null distribution is available in closed form, and says plainly
where an approximation enters. The notation follows “Derivation:
univariate Gest”.

## 1. Why a correction is needed

A point near the boundary of the observation window has no neighbours
outside it. Its observed nearest-neighbour distance is therefore biased
upward, and the raw empirical curve underestimates the true one. All
corrections below trade a small amount of information for reduced bias;
they are different estimators, not different null models. Crucially,
under the mark-permutation null *the locations do not move*, so every
correction’s null distribution is still induced purely by the random
choice of .

Let

denote the distance from location to the boundary of the window , and
let mark the points that are still “at risk” at radius .

## 2. Correction none (raw)

The raw estimator is exactly the statistic of Theorem 1,

so its null mean and variance are given exactly by Theorems 1 and 2.
This is the only fully exact combination in the package: mean **and**
variance, for both and , with no approximation of any kind.

## 3. Correction rs (reduced sample / border)

Only points farther than from the boundary contribute:

Write and .

**Mean of the denominator.** is a simple random sample sum, so

**Mean of the numerator.**

**Second moment of the numerator.** With ,

where, exactly as in Lemma 3,

**Cross moment.**

The reason for the two cases is that, conditional on both and being
marked, a hit for is already guaranteed when the two are within ;
otherwise at least one of ’s own neighbours must be selected. Then

**Combining them.** The package reports the ratio of expectations as the
null mean and a first-order (delta method) variance:

where , and are the denominator moments. The first two moments of are
**exact**; only the way they are combined into a mean and variance for
the ratio is approximate, with bias of order . Radii at which the
population expectation vanishes give no estimate, and are returned as
NA.

## 4. Correction han (Hanisch)

The Hanisch estimator re-weights each event by the area still available
at its radius. With bin edges and the eroded-window areas (the area of
the set of window points with boundary distance greater than ),

where counts the uncensored events whose observed distance falls in bin
. Equivalently, counts marked points with and .

Under the null, conditional on , at least one selected neighbour within
distance has probability

so the *exact* expected bin counts are

The package computes these expectations exactly and then forms the
corrected curve by plug-in, replacing each with . The expected bin
counts are exact; the curve assembled from them is a mean-field
approximation, because for the nonlinear normalisation. No closed-form
variance is available, so csr_var is NA; use the permutation fallback if
a null band for han is needed.

## 5. Correction km (Kaplan-Meier)

The Kaplan-Meier correction treats a point whose nearest neighbour lies
beyond the boundary as right-censored at . The observed variable is
again with event indicator , and the product-limit estimate of the
survival function is converted to a distribution function:

where is the number of uncensored events in bin (so above) and is the
number at risk just before bin .

The expected histogram of is exact as well:

The last term is the censoring mass: when the boundary distance itself
falls in bin , then exactly and no event can occur beyond it. Setting

the reported mean-field curve is . As with han the bin expectations are
exact and the product-limit assembly is a plug-in, so no closed-form
variance is reported.

## 6. What is exact

| correction | null mean | null variance | status |
|----|----|----|----|
| none | exact (Theorem 1 / Theorem 3) | exact (Theorem 2 / Theorem 4) | fully exact |
| rs | ratio of exact expectations | delta method from exact moments | approximate, bias O(1/n) |
| han | plug-in of exact bin expectations | not available | mean-field; permutation for variance |
| km | plug-in of exact bin expectations | not available | mean-field; permutation for variance |

The package records this in the returned object: exact is TRUE only for
correction = “none” with method = “exact”.

## 7. Practical guidance

- Because the window built automatically from a data frame is the
  **convex hull** of the points, many locations sit on or near its
  boundary and can be very small. Edge-corrected analyses are much more
  informative when the user supplies a meaningful window, for example by
  passing a spatstat ppp object with an explicit observation window.
- The reduced-sample estimator is undefined for radii larger than the
  largest ; those rows are returned as NA.
- The Hanisch curve is undefined once the eroded window has zero area,
  and those rows are NA as well.
- If a null band is required for han or km, use permute_gest with the
  same correction.

## 8. Numerical check: rs against exhaustive enumeration

With an explicit unit-square window, every marked subset can be
enumerated for a small pattern and compared with the closed-form
reduced-sample moments.

``` r

suppressPackageStartupMessages(library(permfreeG))
set.seed(3)
N <- 10; n <- 4
d <- data.frame(x = stats::runif(N), y = stats::runif(N))
d$m <- factor(c(rep("A", n), rep("B", N - n)))
pp <- spatstat.geom::ppp(d$x, d$y,
                         window = spatstat.geom::owin(c(0, 1), c(0, 1)),
                         marks = d$m)
radii <- seq(0, 0.3, by = 0.05)

D <- as.matrix(stats::dist(cbind(d$x, d$y)))
bd <- spatstat.geom::bdist.points(pp)
ci <- utils::combn(N, n)
est <- matrix(NA_real_, ncol(ci), length(radii))
for (j in seq_len(ncol(ci))) {
  S <- ci[, j]
  nnd <- apply(D[S, S], 1, function(z) min(z[z > 0]))
  for (k in seq_along(radii)) {
    den <- sum(bd[S] > radii[k])
    est[j, k] <- if (den > 0) {
      sum(nnd <= radii[k] & bd[S] > radii[k]) / den
    } else NA_real_
  }
}
ok <- colSums(!is.na(est)) == ncol(ci)
res <- as.data.frame(exact_gest(pp, radii = radii, marks_i = "A",
                                correction = "rs"))
data.frame(r = radii[ok],
           enumerated_mean = colMeans(est[, ok, drop = FALSE]),
           exact_rs_mean   = res$csr_mean[ok],
           enumerated_var  = apply(est[, ok, drop = FALSE], 2,
                                   function(z) mean((z - mean(z))^2)),
           delta_var       = res$csr_var[ok])
#>      r enumerated_mean exact_rs_mean enumerated_var  delta_var
#> 1 0.00       0.0000000     0.0000000     0.00000000 0.00000000
#> 2 0.05       0.1250000     0.1250000     0.05312500 0.05721726
#> 3 0.10       0.2416667     0.2416667     0.07582341 0.08085813
#> 4 0.15       0.1809524     0.1785714     0.07148904 0.07238737
```

The ratio of expectations reproduces the enumerated mean to within a
fraction of a percent, and the delta-method variance is close to but not
identical with the enumerated variance, as expected.

## 9. Numerical check: han and km against permutation

On a window with mild censoring, the mean-field curves are close to the
permutation average.

``` r

set.seed(99)
N <- 40; n <- 15
d2 <- data.frame(x = stats::runif(N), y = stats::runif(N))
d2$m <- factor(c(rep("A", n), rep("B", N - n)))
pp2 <- spatstat.geom::ppp(d2$x, d2$y,
                          window = spatstat.geom::owin(c(0, 1), c(0, 1)),
                          marks = d2$m)
radii <- seq(0, 0.4, by = 0.05)

out <- lapply(c("han", "km"), function(cor) {
  ex <- as.data.frame(exact_gest(pp2, radii = radii, marks_i = "A",
                                 correction = cor))
  pm <- permute_gest(pp2, radii = radii, marks_i = "A", correction = cor,
                     nsim = 200)
  data.frame(correction = cor, r = radii,
             mean_field = ex$csr_mean, permutation_mean = pm$mean)
})
out <- do.call(rbind, out)
aggregate(abs(out$mean_field - out$permutation_mean),
          by = list(correction = out$correction), FUN = max)
#>   correction          x
#> 1        han 0.01480883
#> 2         km 0.05338997
```

The mean-field approximation is adequate here (differences of a few
percent) but should not be trusted when the marked count is very small
or when most points are heavily censored; in those regimes the
permutation fallback is the right tool.
