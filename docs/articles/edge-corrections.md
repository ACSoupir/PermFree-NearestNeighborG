# Supplement: edge corrections for the permutation-free G function

This vignette supplements the two derivation vignettes. It defines each
edge correction supported by the package, states exactly which part of
its null distribution is available in closed form, and says plainly
where an approximation enters. The notation follows “Derivation:
univariate Gest”.

## 1. Why a correction is needed

A point p near the boundary of the observation window has no neighbours
outside it. Its observed nearest-neighbour distance is therefore biased
upward, and the raw empirical G curve underestimates the true one. All
corrections below trade a small amount of information for reduced bias;
they are different estimators, not different null models. Crucially,
under the mark-permutation null *the locations do not move*, so every
correction’s null distribution is still induced purely by the random
choice of S.

Let

b_p = d(p, W^c)

denote the distance from location p to the boundary of the window W, and
let u_p(r) = \mathbf{1}\\b_p \> r\\ mark the points that are still “at
risk” at radius r.

## 2. Correction none (raw)

The raw estimator is exactly the statistic of Theorem 1,

\hat{G}\_{\mathrm{none}}(r) = G_S(r) = \frac{1}{n}\sum\_{p\in S} X_p,

so its null mean and variance are given exactly by Theorems 1 and 2.
This is the only fully exact combination in the package: mean **and**
variance, for both G and G\_{ij}, with no approximation of any kind.

## 3. Correction rs (reduced sample / border)

Only points farther than r from the boundary contribute:

\hat{G}\_{\mathrm{rs}}(r) = \frac{Z(r)}{Y(r)}, \qquad Y(r) =
\sum\_{p=1}^{N} u_p(r)\\\mathbf{1}\\p \in S\\, \qquad Z(r) =
\sum\_{p=1}^{N} u_p(r)\\ X_p .

Write \bar{u}(r) = N^{-1}\sum_p u_p(r) and a_p = 1 - T\_{N-1,n-1}(m_p).

**Mean of the denominator.** Y is a simple random sample sum, so

\mathbb{E}\[Y\] = n\\\bar{u}, \qquad \mathrm{Var}(Y) =
\frac{n(N-n)}{N-1}\\\bar{u}\\(1-\bar{u}).

**Mean of the numerator.**

\mathbb{E}\[Z\] = \frac{n}{N}\sum\_{p=1}^{N} u_p\\ a_p .

**Second moment of the numerator.** With c_2 = n(n-1)/\[N(N-1)\],

\mathbb{E}\[Z^2\] = \sum_p u_p \frac{n}{N} a_p + 2\sum\_{p\<q} u_p u_q\\
c_2\\ g^{\mathrm{both}}\_{pq},

where, exactly as in Lemma 3,

g^{\mathrm{both}}\_{pq} = \begin{cases} 1, & \\x_p-x_q\\ \le r,\\ 1 -
T\_{N-2,n-2}(m_p) - T\_{N-2,n-2}(m_q) + T\_{N-2,n-2}(m_p+m_q-w\_{pq}), &
\text{otherwise}. \end{cases}

**Cross moment.**

\mathbb{E}\[ZY\] = \sum_p u_p\frac{n}{N}a_p + 2\sum\_{p\<q}u_pu_q\\c_2\\
g^{\mathrm{one}}\_{pq}, \qquad g^{\mathrm{one}}\_{pq} = \begin{cases} 1,
& \\x_p-x_q\\ \le r,\\ 1 - T\_{N-2,n-2}(m_p), & \text{otherwise}.
\end{cases}

The reason for the two cases is that, conditional on both p and q being
marked, a hit for p is already guaranteed when the two are within r;
otherwise at least one of p’s own neighbours must be selected. Then

\mathrm{Cov}(Z,Y) = \mathbb{E}\[ZY\] - \mathbb{E}\[Z\]\mathbb{E}\[Y\].

**Combining them.** The package reports the ratio of expectations as the
null mean and a first-order (delta method) variance:

\mathbb{E}\\\left\[\hat{G}\_{\mathrm{rs}}\right\] \\\approx\\
\frac{\mu_Z}{\mu_Y},

\mathrm{Var}\\\left(\hat{G}\_{\mathrm{rs}}\right) \\\approx\\
\frac{V_Z}{\mu_Y^{2}} - \frac{2\mu_Z\\\mathrm{Cov}(Z,Y)}{\mu_Y^{3}} +
\frac{\mu_Z^{2} V_Y}{\mu_Y^{4}},

where \mu_Z = \mathbb{E}\[Z\], V_Z = \mathrm{Var}(Z) and \mu_Y, V_Y are
the denominator moments. The first two moments of (Z, Y) are **exact**;
only the way they are combined into a mean and variance for the ratio is
approximate, with bias of order 1/n. Radii at which the population
expectation \mu_Y vanishes give no estimate, and are returned as NA.

## 4. Correction han (Hanisch)

The Hanisch estimator re-weights each event by the area still available
at its radius. With bin edges 0 = r_0 \< r_1 \< \dots and the
eroded-window areas a(r) = \|W \ominus r\| (the area of the set of
window points with boundary distance greater than r),

\hat{G}\_{\mathrm{han}}(r_k) = \frac{\sum\_{j\le k} h_j / a(r_j)}
{\max_k \sum\_{j\le k} h_j/a(r_j)},

where h_j counts the uncensored events whose observed distance o_p =
\min(\mathrm{nnd}\_p, b_p) falls in bin j. Equivalently, h_j counts
marked points with \mathrm{nnd}\_p \in (r\_{j-1}, r_j\] and
\mathrm{nnd}\_p \le b_p.

Under the null, conditional on p \in S, at least one selected neighbour
within distance t has probability

F_p(t) = 1 - T\_{N-1,n-1}\big(m_p(t)\big), \qquad m_p(t) = \|\\q \neq p
: \\x_p-x_q\\ \le t\\\|,

so the *exact* expected bin counts are

\mathbb{E}\[h_j\] = \frac{n}{N}\sum\_{p=1}^{N} \Big\[
F_p\\\left(\min(r_j, b_p)\right) - F_p\\\left(\min(r\_{j-1}, b_p)\right)
\Big\].

The package computes these expectations exactly and then forms the
corrected curve by plug-in, replacing each h_j with \mathbb{E}\[h_j\].
The expected bin counts are exact; the curve assembled from them is a
mean-field approximation, because \mathbb{E}\[f(H)\] \neq f(\mathbb{E}H)
for the nonlinear normalisation. No closed-form variance is available,
so csr_var is NA; use the permutation fallback if a null band for han is
needed.

## 5. Correction km (Kaplan-Meier)

The Kaplan-Meier correction treats a point whose nearest neighbour lies
beyond the boundary as right-censored at b_p. The observed variable is
again o_p = \min(\mathrm{nnd}\_p, b_p) with event indicator \delta_p =
\mathbf{1}\\\mathrm{nnd}\_p \le b_p\\, and the product-limit estimate of
the survival function is converted to a distribution function:

\hat{G}\_{\mathrm{km}}(r_k) = 1 - \prod\_{j\le k} \left(1 -
\frac{\mathrm{nco}\_j}{d_j}\right),

where \mathrm{nco}\_j is the number of uncensored events in bin j (so
\mathbb{E}\[\mathrm{nco}\_j\] = \mathbb{E}\[h_j\] above) and d_j is the
number at risk just before bin j.

The expected histogram of o_p is exact as well:

\mathbb{E}\[\mathrm{obs}\_j\] = \frac{n}{N}\sum\_{p=1}^{N} \Big\[
F_p(\min(r_j,b_p)) - F_p(\min(r\_{j-1},b_p)) + \mathbf{1}\\r\_{j-1} \<
b_p \le r_j\\\big(1 - F_p(b_p)\big) \Big\].

The last term is the censoring mass: when the boundary distance itself
falls in bin j, then o_p = b_p exactly and no event can occur beyond it.
Setting

d_k = \sum\_{j \ge k} \mathbb{E}\[\mathrm{obs}\_j\], \qquad s_k = 1 -
\frac{\mathbb{E}\[\mathrm{nco}\_k\]}{d_k},

the reported mean-field curve is 1 - \prod\_{j\le k} s_j. As with han
the bin expectations are exact and the product-limit assembly is a
plug-in, so no closed-form variance is reported.

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
  boundary and b_p can be very small. Edge-corrected analyses are much
  more informative when the user supplies a meaningful window, for
  example by passing a spatstat ppp object with an explicit observation
  window.
- The reduced-sample estimator is undefined for radii larger than the
  largest b_p; those rows are returned as NA.
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
