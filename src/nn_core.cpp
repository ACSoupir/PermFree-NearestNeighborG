#include <Rcpp.h>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Core utilities
// ---------------------------------------------------------------------------

namespace {

inline int popcount64(std::uint64_t x) { return __builtin_popcountll(x); }

// Index of the pair (i, j) with i < j in a packed lower triangle.
inline std::size_t pidx(int N, int i, int j) {
  return static_cast<std::size_t>(i) * N -
         (static_cast<std::size_t>(i) * (i + 1)) / 2 +
         static_cast<std::size_t>(j - i - 1);
}

// T(k) = choose(B - k, s) / choose(B, s), 0 when B - k < s.
std::vector<double> make_table(int B, int s) {
  std::vector<double> T(static_cast<std::size_t>(B) + 1, 0.0);
  if (s < 0 || B - s < 0) return T;
  double log_den = std::lgamma(static_cast<double>(B) + 1.0) -
                   std::lgamma(static_cast<double>(s) + 1.0) -
                   std::lgamma(static_cast<double>(B - s) + 1.0);
  for (int k = 0; k <= B; ++k) {
    if (B - k < s) { T[k] = 0.0; continue; }
    double log_num = std::lgamma(static_cast<double>(B - k) + 1.0) -
                     std::lgamma(static_cast<double>(s) + 1.0) -
                     std::lgamma(static_cast<double>(B - k - s) + 1.0);
    T[k] = std::exp(log_num - log_den);
  }
  return T;
}

inline double tget(const std::vector<double>& T, int k) {
  if (k < 0 || static_cast<std::size_t>(k) >= T.size()) return 0.0;
  return T[static_cast<std::size_t>(k)];
}

// Packed lower-triangle of squared distances.
struct PairDist {
  int N;
  std::vector<double> d2;

  PairDist(const NumericVector& x, const NumericVector& y) : N(x.size()) {
    d2.resize(static_cast<std::size_t>(N) * static_cast<std::size_t>(std::max(N - 1, 0)) / 2);
    for (int i = 0; i < N - 1; ++i) {
      double xi = x[i], yi = y[i];
      for (int j = i + 1; j < N; ++j) {
        double dx = x[j] - xi, dy = y[j] - yi;
        d2[pidx(N, i, j)] = dx * dx + dy * dy;
      }
    }
  }
};

// Neighbour bitsets and neighbour counts at squared radius r2.
inline void build_bits(const PairDist& P, double r2, int W,
                       std::vector<std::uint64_t>& B, std::vector<int>& m) {
  const int N = P.N;
  std::fill(B.begin(), B.end(), static_cast<std::uint64_t>(0));
  std::fill(m.begin(), m.end(), 0);
  for (int i = 0; i < N - 1; ++i) {
    std::uint64_t* bi = &B[static_cast<std::size_t>(i) * W];
    for (int j = i + 1; j < N; ++j) {
      if (P.d2[pidx(N, i, j)] <= r2) {
        bi[j >> 6] |= (1ULL << (j & 63));
        B[static_cast<std::size_t>(j) * W + (i >> 6)] |= (1ULL << (i & 63));
        m[i]++;
        m[j]++;
      }
    }
  }
}

struct Moments {
  std::vector<double> mean;
  std::vector<double> var;
};

// Exact CSR moments for the permutation null.
//
// cross = false : subset of size n_i drawn uniformly from all N locations
//                 (univariate Gest).
// cross = true  : disjoint label permutation, type-i set size n_i and
//                 type-j set size n_j (bivariate Gcross).
Moments moments_core(const PairDist& P, const NumericVector& radii,
                     int n_i, bool cross, int n_j, bool covariance) {
  const int N = P.N;
  const int Rr = radii.size();
  Moments out;
  out.mean.assign(Rr, NA_REAL);
  out.var.assign(Rr, NA_REAL);
  if (N < 2 || n_i < 1) return out;

  const int W = (N + 63) / 64;
  std::vector<std::uint64_t> B(static_cast<std::size_t>(N) * W);
  std::vector<int> m(N);

  const int s_eff = cross ? n_j : (n_i - 1);
  const int t_eff = cross ? n_j : (n_i - 2);
  std::vector<double> T1 = make_table(N - 1, s_eff);
  std::vector<double> T2;
  if (covariance) T2 = make_table(N - 2, t_eff);

  const double c2 = static_cast<double>(n_i) * (n_i - 1) /
                    (static_cast<double>(N) * (N - 1));
  const double n_over_N = static_cast<double>(n_i) / N;
  std::vector<double> EX(N);

  for (int k = 0; k < Rr; ++k) {
    double r2 = static_cast<double>(radii[k]) * static_cast<double>(radii[k]);
    build_bits(P, r2, W, B, m);

    double mu = 0.0;
    for (int p = 0; p < N; ++p) {
      double a_p = 1.0 - tget(T1, m[p]);
      mu += a_p;
      EX[p] = n_over_N * a_p;
    }
    out.mean[k] = mu / N;

    double vv = 0.0;
    for (int p = 0; p < N; ++p) vv += EX[p] * (1.0 - EX[p]);

    if (covariance) {
      for (int p = 0; p < N - 1; ++p) {
        const std::uint64_t* bp = &B[static_cast<std::size_t>(p) * W];
        const int mp = m[p];
        for (int q = p + 1; q < N; ++q) {
          bool nb = ((bp[q >> 6] >> (q & 63)) & 1ULL) != 0ULL;
          double Eb;
          if (cross) {
            const int mps = mp - (nb ? 1 : 0);
            const int mqs = m[q] - (nb ? 1 : 0);
            const std::uint64_t* bq = &B[static_cast<std::size_t>(q) * W];
            int wc = 0;
            for (int wi = 0; wi < W; ++wi) wc += popcount64(bp[wi] & bq[wi]);
            const int un = mps + mqs - wc;
            Eb = c2 * (1.0 - tget(T2, mps) - tget(T2, mqs) + tget(T2, un));
          } else if (nb) {
            Eb = c2;
          } else {
            const std::uint64_t* bq = &B[static_cast<std::size_t>(q) * W];
            int wc = 0;
            for (int wi = 0; wi < W; ++wi) wc += popcount64(bp[wi] & bq[wi]);
            const int un = mp + m[q] - wc;
            Eb = c2 * (1.0 - tget(T2, mp) - tget(T2, m[q]) + tget(T2, un));
          }
          vv += 2.0 * (Eb - EX[p] * EX[q]);
        }
      }
    }

    double var = vv / (static_cast<double>(n_i) * n_i);
    out.var[k] = var > 0.0 ? var : 0.0;
  }
  return out;
}

} // namespace

// ---------------------------------------------------------------------------
// Exported functions
// ---------------------------------------------------------------------------

//' Neighbour counts within each radius (<= r)
//'
//' @param x,y Point coordinates.
//' @param radii Radii (must be increasing for the binary-search fast path; any
//'   order is accepted).
//' @return An \code{N x length(radii)} integer matrix.
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix nn_counts_cpp(NumericVector x, NumericVector y, NumericVector radii) {
  const int N = x.size();
  const int Rr = radii.size();
  IntegerMatrix out(N, Rr);
  for (int p = 0; p < N; ++p) {
    std::vector<double> d2(N);
    for (int q = 0; q < N; ++q) {
      double dx = x[q] - x[p], dy = y[q] - y[p];
      d2[q] = dx * dx + dy * dy;
    }
    for (int k = 0; k < Rr; ++k) {
      double t2 = static_cast<double>(radii[k]) * static_cast<double>(radii[k]);
      int c = 0;
      for (int q = 0; q < N; ++q) {
        if (q != p && d2[q] <= t2) c++;
      }
      out(p, k) = c;
    }
  }
  return out;
}

//' Neighbour counts within min(radius, boundary distance)
//'
//' Used by the Hanisch and Kaplan-Meier moment calculations: an event is only
//' observed when its nearest-neighbour distance does not exceed the point's
//' distance to the window boundary.
//'
//' @param x,y Point coordinates.
//' @param bdist Distance of each point to the window boundary.
//' @param radii Radii.
//' @return An \code{N x length(radii)} integer matrix with entry
//'   \eqn{\#\{q \ne p : d_{pq} \le \min(r_k, b_p)\}}.
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix censored_neighbor_counts_cpp(NumericVector x, NumericVector y,
                                           NumericVector bdist,
                                           NumericVector radii) {
  const int N = x.size();
  const int Rr = radii.size();
  IntegerMatrix out(N, Rr);
  for (int p = 0; p < N; ++p) {
    std::vector<double> d2(N);
    for (int q = 0; q < N; ++q) {
      double dx = x[q] - x[p], dy = y[q] - y[p];
      d2[q] = dx * dx + dy * dy;
    }
    for (int k = 0; k < Rr; ++k) {
      double t = std::min(static_cast<double>(radii[k]),
                          static_cast<double>(bdist[p]));
      double t2 = t * t;
      int c = 0;
      for (int q = 0; q < N; ++q) {
        if (q != p && d2[q] <= t2) c++;
      }
      out(p, k) = c;
    }
  }
  return out;
}

//' Exact univariate CSR mean and variance
//'
//' @param x,y Point coordinates of the full pattern.
//' @param radii Radii at which to evaluate the moments.
//' @param n_subset Number of marked points selected in each random subset.
//' @param covariance Include the off-diagonal covariance terms (exact variance)
//'   or only the marginal term (independence approximation).
//' @return A list with numeric vectors \code{mean} and \code{var}.
//' @keywords internal
// [[Rcpp::export]]
List csr_uni_moments_cpp(NumericVector x, NumericVector y,
                         NumericVector radii, int n_subset,
                         bool covariance = true) {
  PairDist P(x, y);
  Moments M = moments_core(P, radii, n_subset, false, -1, covariance);
  return List::create(_["mean"] = wrap(M.mean), _["var"] = wrap(M.var));
}

//' Exact bivariate Gcross CSR mean and variance
//'
//' Disjoint label-permutation null: type-i set of size \code{n_i}, type-j set
//' of size \code{n_j} drawn without replacement from the remaining locations.
//'
//' @param x,y Point coordinates of the full pattern.
//' @param radii Radii at which to evaluate the moments.
//' @param n_i Number of type-i points in each permutation.
//' @param n_j Number of type-j points in each permutation.
//' @param covariance Include the off-diagonal covariance terms (exact variance)
//'   or only the marginal term.
//' @return A list with numeric vectors \code{mean} and \code{var}.
//' @keywords internal
// [[Rcpp::export]]
List csr_cross_moments_cpp(NumericVector x, NumericVector y,
                           NumericVector radii, int n_i, int n_j,
                           bool covariance = true) {
  PairDist P(x, y);
  Moments M = moments_core(P, radii, n_i, true, n_j, covariance);
  return List::create(_["mean"] = wrap(M.mean), _["var"] = wrap(M.var));
}

//' Reduced-sample (border) CSR moments
//'
//' Exact first and second moments of the numerator Z and denominator Y of the
//' reduced-sample estimator \eqn{\hat G_{rs}(r) = Z/Y}, combined into the
//' ratio-of-expectations mean and a first-order (delta method) variance.
//'
//' @param x,y Point coordinates of the full pattern.
//' @param bdist Distance from each point to the window boundary.
//' @param radii Radii at which to evaluate the moments.
//' @param n_i Size of the randomly selected subset (type-i set for cross).
//' @param n_j For bivariate Gcross the size of the type-j set; use \code{-1}
//'   for univariate Gest.
//' @return A list with \code{mean}, \code{var} and the raw components
//'   \code{z_mean}, \code{z_var}, \code{y_mean}, \code{y_var}, \code{cov}.
//' @keywords internal
// [[Rcpp::export]]
List rs_moments_cpp(NumericVector x, NumericVector y, NumericVector bdist,
                    NumericVector radii, int n_i, int n_j = -1) {
  const bool cross = (n_j >= 0);
  PairDist P(x, y);
  const int N = P.N;
  const int Rr = radii.size();
  NumericVector mean(Rr, NA_REAL), var(Rr, NA_REAL);
  NumericVector z_mean(Rr, NA_REAL), z_var(Rr, NA_REAL);
  NumericVector y_mean(Rr, NA_REAL), y_var(Rr, NA_REAL), cov(Rr, NA_REAL);
  if (N < 2 || n_i < 1) {
    return List::create(_["mean"] = mean, _["var"] = var,
                        _["z_mean"] = z_mean, _["z_var"] = z_var,
                        _["y_mean"] = y_mean, _["y_var"] = y_var, _["cov"] = cov);
  }

  const int W = (N + 63) / 64;
  std::vector<std::uint64_t> B(static_cast<std::size_t>(N) * W);
  std::vector<int> m(N);

  const int s_eff = cross ? n_j : (n_i - 1);
  const int t_eff = cross ? n_j : (n_i - 2);
  std::vector<double> T1 = make_table(N - 1, s_eff);
  std::vector<double> T2 = make_table(N - 2, t_eff);

  const double c2 = static_cast<double>(n_i) * (n_i - 1) /
                    (static_cast<double>(N) * (N - 1));
  const double n_over_N = static_cast<double>(n_i) / N;
  std::vector<double> u(N), a_p(N);

  for (int k = 0; k < Rr; ++k) {
    double r = radii[k];
    build_bits(P, r * r, W, B, m);

    double su = 0.0;
    for (int p = 0; p < N; ++p) {
      u[p] = (bdist[p] > r) ? 1.0 : 0.0;
      su += u[p];
      a_p[p] = 1.0 - tget(T1, m[p]);
    }
    const double EY = n_over_N * su;
    if (EY <= 0.0) continue;

    double EZ = 0.0;
    for (int p = 0; p < N; ++p) EZ += u[p] * n_over_N * a_p[p];

    double EZ2 = 0.0, EZY = 0.0;
    for (int p = 0; p < N; ++p) {
      double diag = u[p] * n_over_N * a_p[p];
      EZ2 += diag;
      EZY += diag;
    }
    for (int p = 0; p < N - 1; ++p) {
      if (u[p] == 0.0) continue;
      const std::uint64_t* bp = &B[static_cast<std::size_t>(p) * W];
      for (int q = p + 1; q < N; ++q) {
        if (u[q] == 0.0) continue;
        bool nb = ((bp[q >> 6] >> (q & 63)) & 1ULL) != 0ULL;
        const int mp = m[p], mq = m[q];
        double p_both, p_p;
        if (cross) {
          const int mps = mp - (nb ? 1 : 0);
          const std::uint64_t* bq = &B[static_cast<std::size_t>(q) * W];
          int wc = 0;
          for (int wi = 0; wi < W; ++wi) wc += popcount64(bp[wi] & bq[wi]);
          p_both = c2 * (1.0 - tget(T2, mps) -
                         tget(T2, mq - (nb ? 1 : 0)) +
                         tget(T2, mps + (mq - (nb ? 1 : 0)) - wc));
          p_p = c2 * (1.0 - tget(T2, mps));
        } else if (nb) {
          p_both = c2;
          p_p = c2;
        } else {
          const std::uint64_t* bq = &B[static_cast<std::size_t>(q) * W];
          int wc = 0;
          for (int wi = 0; wi < W; ++wi) wc += popcount64(bp[wi] & bq[wi]);
          p_both = c2 * (1.0 - tget(T2, mp) - tget(T2, mq) +
                         tget(T2, mp + mq - wc));
          p_p = c2 * (1.0 - tget(T2, mp));
        }
        const double ww = 2.0 * u[p] * u[q];
        EZ2 += ww * p_both;
        EZY += ww * p_p;
      }
    }

    double varZ = EZ2 - EZ * EZ; if (varZ < 0.0) varZ = 0.0;
    const double ubar = su / N;
    double vy = (N > 1)
      ? static_cast<double>(n_i) * (N - n_i) / (N - 1) * ubar * (1.0 - ubar)
      : 0.0;
    const double cZY = EZY - EZ * EY;

    z_mean[k] = EZ; z_var[k] = varZ;
    y_mean[k] = EY; y_var[k] = vy; cov[k] = cZY;
    mean[k] = EZ / EY;
    double vratio = varZ / (EY * EY)
                  - 2.0 * EZ * cZY / (EY * EY * EY)
                  + (EZ * EZ) * vy / (EY * EY * EY * EY);
    var[k] = vratio > 0.0 ? vratio : 0.0;
  }

  return List::create(_["mean"] = mean, _["var"] = var,
                      _["z_mean"] = z_mean, _["z_var"] = z_var,
                      _["y_mean"] = y_mean, _["y_var"] = y_var, _["cov"] = cov);
}

//' Number of neighbours closer than the window boundary
//'
//' For each point p, the number of other points q whose distance to p is at
//' most the boundary distance b_p.  Used by the
//' Kaplan-Meier / Hanisch mean-field moment calculations.
//'
//' @param x,y Point coordinates.
//' @param bdist Distance of each point to the window boundary.
//' @return Integer vector with one entry per point.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector counts_within_bdist_cpp(NumericVector x, NumericVector y,
                                      NumericVector bdist) {
  const int N = x.size();
  IntegerVector out(N);
  for (int p = 0; p < N; ++p) {
    double b2 = static_cast<double>(bdist[p]) * static_cast<double>(bdist[p]);
    int c = 0;
    for (int q = 0; q < N; ++q) {
      if (q == p) continue;
      double dx = x[q] - x[p], dy = y[q] - y[p];
      if (dx * dx + dy * dy <= b2) c++;
    }
    out[p] = c;
  }
  return out;
}
