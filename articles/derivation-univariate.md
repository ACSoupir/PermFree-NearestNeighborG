# Derivation: exact mean and variance of the univariate nearest-neighbour G function

This vignette proves, without reference to any code, the closed-form
mean and variance of the univariate nearest-neighbour function under the
mark-permutation null. The bivariate case is proved in the companion
vignette “Derivation: bivariate Gcross”, and the edge corrections are
derived as a supplement in “Edge corrections”.

## 1. The null model and the statistic

Let be the observed spatial locations. These are held **fixed**. A mark
(for example a positive cell phenotype) is assigned to exactly of them,
and the null hypothesis of complete spatial randomness (CSR) for the
marks says that every assignment is equally likely. Equivalently, the
marked set

is a uniformly random -subset of the locations.

For write for a point and a set , and define the sample
nearest-neighbour distance of within as . The univariate statistic is
the empirical distribution function of those distances over the marked
set,

Throughout we use the closed convention “”, which is exactly what
spatstat.explore::Gest() computes; the difference from a strict
inequality is a null set for continuous locations.

## 2. Indicator decomposition

The key step replaces the average over by a sum of point indicators. For
each location put

and

Then

Indeed, the term for is one precisely when belongs to the marked set and
some other selected point lies within distance ; dividing by turns the
count into the proportion of marked points that have a neighbour.

## 3. The avoidance table

All moments below are built from one combinatorial quantity. For
integers and , define

is the probability that a uniformly drawn -subset of a set with elements
avoids prescribed elements. The package precomputes it once per analysis
from logarithms of the gamma function, which is why no binomial
coefficients are ever formed explicitly.

**Lemma 1 (avoidance).** Let be a uniform -subset of a set with , and
let with . Then

*Proof.* Choosing uniformly is equivalent to choosing its complement , a
uniform -subset, and holds exactly when that complement contains all of
. The number of -subsets containing a fixed -set is ; dividing by gives
the claim. Both sides are zero when .

## 4. Exact mean

**Theorem 1 (mean).** With the notation of Section 2,

*Proof.* By linearity and the decomposition of Section 2,

Fix . Since , the multiplication rule gives

Conditional on , the remaining members of form a uniform -subset of the
locations in , and is a fixed subset of that ground set with elements.
Lemma 1 therefore gives

Substituting and using completes the proof.

Two immediate consequences are worth recording.

- The mean is **sample-specific**: it depends on the observed locations
  through the neighbour counts alone.
- It never depends on which locations happen to carry the mark, only on
  ; relabelling the marks leaves unchanged.

## 5. Exact second moments

For the variance we need both marginal and joint probabilities.

**Lemma 2 (pair inclusion).** For distinct locations ,

*Proof.* There are equally likely marked sets; those containing both and
number . The ratio is .

**Lemma 3 (joint emptiness).** Let and suppose , so that and . Write .
Then

*Proof.* Conditional on , the other marked locations form a uniform
-subset of the remaining locations. Both and are fixed subsets of that
ground set, of sizes and , with intersection of size . The event that
meets neither set is the complement of “meets ” or “meets ”, and the
union has size . Applying Lemma 1 to each of the three sets gives the
stated expression.

Combining Lemmas 2 and 3,

The first case is immediate: if itself lies within distance of , then
selecting both already guarantees that each has a neighbour, so the
conditional probability of two hits is one.

**Theorem 2 (variance).**

where and is given above.

*Proof.* Expand the square of the sum:

because the indicators are Bernoulli, so , and covariance is bilinear.

The result is **exact** for every finite sample: no asymptotic
approximation, no independence assumption and no simulation are used.
The only numerical input is the avoidance table, computed in log space.

## 6. Computational cost

Computing all costs one pass over the pairs per radius. The double sum
has one term per pair as well, but each needs the intersection size .
Storing each as a packed bit set of machine words and using a population
count makes each intersection rather than , so the total work is

for radii. The original implementation instead called a binomial
coefficient routine inside the inner loop, which is what made the exact
variance prohibitive for large .

If only a rough scale is needed, setting covariance = FALSE keeps the
diagonal term and drops ; this is the “independent points”
approximation. It costs after the neighbour counts are known.

## 7. Numerical illustration

The following small example checks Theorem 1 and Theorem 2 against
*exhaustive enumeration* of all marked sets.

``` r

suppressPackageStartupMessages(library(permfreeG))
set.seed(42)
N <- 10; n <- 4
x <- stats::runif(N); y <- stats::runif(N)
radii <- seq(0, 1.2, by = 0.15)

# exhaustive enumeration over all C(N, n) subsets
D <- as.matrix(stats::dist(cbind(x, y)))
ci <- utils::combn(N, n)
Gsim <- matrix(NA_real_, ncol(ci), length(radii))
for (j in seq_len(ncol(ci))) {
  S <- ci[, j]
  DS <- D[S, S]; diag(DS) <- Inf
  Gsim[j, ] <- vapply(radii, function(k) mean(rowSums(DS <= k) > 0), numeric(1))
}

res <- as.data.frame(exact_gest(data.frame(x = x, y = y,
                                           m = factor(c(rep("A", n),
                                                        rep("B", N - n)))),
                                radii = radii, mark_col = "m", marks_i = "A"))

max_mean_err <- max(abs(res$csr_mean - colMeans(Gsim)))
max_var_err  <- max(abs(res$csr_var -
                        apply(Gsim, 2, function(z) mean((z - mean(z))^2))))
c(max_mean_error = max_mean_err, max_variance_error = max_var_err)
#>     max_mean_error max_variance_error 
#>       4.440892e-16       1.822544e-16
```

Both errors are at the level of floating-point round-off, confirming
that Theorems 1 and 2 reproduce the permutation distribution exactly.

``` r

plot(radii, colMeans(Gsim), type = "l", lwd = 3, col = "grey80",
     xlab = "r", ylab = "G(r)", main = "Exact null vs exhaustive enumeration")
lines(res$r, res$csr_mean, col = "red", lty = 2)
sd <- sqrt(res$csr_var)
lines(res$r, res$csr_mean - sd, col = "red", lty = 3)
lines(res$r, res$csr_mean + sd, col = "red", lty = 3)
legend("bottomright", c("enumeration mean", "exact mean", "+/- 1 sd"),
       col = c("grey80", "red", "red"), lty = c(1, 2, 3), bty = "n")
```

![Exact CSR mean and standard deviation band overlaid on the
exhaustive-enumeration
mean.](derivation-univariate_files/figure-html/plot-1.png)
