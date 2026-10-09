#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_neighbor_counts(NumericVector x, NumericVector y, NumericVector radii) {
  int n = x.size();
  int nr = radii.size();
  NumericMatrix counts(n, nr);
  
  // Precompute squared radii to avoid repeated sqrt computations.
  NumericVector radii2(nr);
  for(int k = 0; k < nr; k++){
    radii2[k] = radii[k] * radii[k];
  }
  
  // Loop over all points
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i == j) continue;  // Skip self-comparison
      double dx = x[j] - x[i];
      double dy = y[j] - y[i];
      double d2 = dx * dx + dy * dy;
      // For each radius, check if the point j is within the circle centered at point i
      for(int k = 0; k < nr; k++){
        if(d2 <= radii2[k]){
          counts(i, k) += 1;
        }
      }
    }
  }
  return counts;
}
