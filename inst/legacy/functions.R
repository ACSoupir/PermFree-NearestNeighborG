
# Functions ---------------------------------------------------------------
message("Compiling C++")
Rcpp::sourceCpp("compute_neighbor_counts.cpp")
Rcpp::sourceCpp("compute_neighbor_counts_opt.cpp")
#faster matrix mulipliation
Rcpp::sourceCpp("mat_mult_eigen.cpp")


message("R functions")
permutations = function(data, radii, mark_col, mark, perms){
  win = spatstat.geom::convexhull.xy(data[,1:2])
  pp = spatstat.geom::ppp(x = data$x,
                          y = data$y,
                          window = win,
                          marks = as.factor(data[[mark_col]]))
  G_obs = spatstat.explore::Gest(subset(pp, marks == mark), r = radii, rmax = max(radii), correction = "none") %>%
    as.data.frame() %>%
    dplyr::rename("Theoretical CSR" = 2, "Observed G" = 3)
  
  pcells = sum(data[[mark_col]] == mark)
  perm = lapply(1:perms, function(x){
    spatstat.explore::Gest(pp[sample(1:nrow(data), pcells, replace = F)], r =radii, rmax = max(radii), correction = "none") %>%
      as.data.frame() %>%
      dplyr::mutate(Group = "Permuted",
                    rep = x)
  }) %>%
    do.call(dplyr::bind_rows, .)%>%
    dplyr::rename("Theoretical CSR" = 2, "Empirical CSR" = 3)
  res = perm %>%
    dplyr::select(1,3,4,5) %>%
    dplyr::full_join(G_obs, ., dplyr::join_by(r))
  return(res)
}

#combinatorial calculations in log space for large N samples
log_comb = function(n, k) {
  lgamma(n + 1) - (lgamma(k + 1) + lgamma(n - k + 1))
}

binom_prob = function(n, k, x) {
  exp(log_comb(n - x, k) - log_comb(n, k))
}

#means
exact_csr = function(data, radii, mark_col, mark){
  win = spatstat.geom::convexhull.xy(data[,1:2])
  pp = spatstat.geom::ppp(x = data$x,
                          y = data$y,
                          window = win,
                          marks = as.factor(data[[mark_col]]))
  G_obs = spatstat.explore::Gest(subset(pp, marks == mark), r = radii, rmax = max(radii), correction = "none") %>%
    as.data.frame() %>%
    dplyr::rename("Theoretical CSR" = 2, "Observed G" = 3)
  
  pcells = sum(data[[mark_col]] == mark)
  summation = sapply(1:nrow(data), function(i){
    dx = data$x[-i] - data$x[i]
    dy = data$y[-i] - data$y[i]
    distances = sqrt(dx^2 + dy^2)
    ni_r = sapply(radii, function(x) sum(distances <= x))
    binom_prob(nrow(data)-1, pcells-1, ni_r)
  })
  
  summation2 = rowSums(summation)
  res = data.frame(r = radii, exact_csr = 1 - 1/nrow(data) * summation2)  %>%
    dplyr::mutate(Group = "Exact") %>%
    dplyr::full_join(G_obs, ., by = dplyr::join_by(r)) %>%
    dplyr::rename("Empirical CSR" = 4)
}

exact_csr_rcpp = function(data, radii, mark_col, mark) {
  win = spatstat.geom::convexhull.xy(data[, 1:2])
  pp = spatstat.geom::ppp(x = data$x,
                           y = data$y,
                           window = win,
                           marks = as.factor(data[[mark_col]]))
  G_obs = spatstat.explore::Gest(subset(pp, marks == mark), r = radii, rmax = max(radii), correction = "none") %>%
    as.data.frame() %>%
    dplyr::rename("Theoretical CSR" = 2, "Observed G" = 3)
  
  pcells = sum(data[[mark_col]] == mark)
  n_points = nrow(data)
  
  counts = compute_neighbor_counts(data$x, data$y, radii)
  binom_vals = binom_prob(n_points-1, pcells-1, counts)
  summation2 = colSums(binom_vals)
  
  exact_csr_values = 1 - (1 / n_points) * summation2
  
  res = data.frame(r = radii, exact_csr = exact_csr_values) %>%
    dplyr::mutate(Group = "Exact") %>%
    dplyr::full_join(G_obs, ., by = dplyr::join_by(r)) %>%
    dplyr::rename("Empirical CSR" = exact_csr)
  
  return(res)
}

exact_csr_rcpp2 = function(data, radii, mark_col, mark) {
  win = spatstat.geom::convexhull.xy(data[, 1:2])
  pp = spatstat.geom::ppp(x = data$x,
                           y = data$y,
                           window = win,
                           marks = as.factor(data[[mark_col]]))
  G_obs = spatstat.explore::Gest(subset(pp, marks == mark), r = radii, rmax = max(radii), correction = "none") %>%
    as.data.frame() %>%
    dplyr::rename("Theoretical CSR" = 2, "Observed G" = 3)
  
  pcells = sum(data[[mark_col]] == mark)
  n_points = nrow(data)
  
  counts = compute_neighbor_counts_opt(data$x, data$y, radii)
  binom_vals = binom_prob(n_points-1, pcells-1, counts)
  summation2 = colSums(binom_vals)
  
  exact_csr_values = 1 - (1 / n_points) * summation2
  
  res = data.frame(r = radii, exact_csr = exact_csr_values) %>%
    dplyr::mutate(Group = "Exact") %>%
    dplyr::full_join(G_obs, ., by = dplyr::join_by(r)) %>%
    dplyr::rename("Empirical CSR" = exact_csr)
  
  return(res)
}

#variance

exact_G_variance = function(df, mark_col, mark, radii) {
  if(nrow(df) > 250) message("Large number of points - this may take a while...")
  N  = nrow(df)
  n  = sum(df[[mark_col]] == mark)
  if (n < 2) stop("Need at least two marked points.")
  ell = n - 2
  D   = as.matrix(dist(df[, c("x", "y")]))
  
  res = numeric(length(radii))
  #names(res) = format(radii)
  
  for (i in seq_along(radii)) {
    r = radii[i]
    neigh  = (D < r) & (D > 0)
    n_p    = rowSums(neigh)
    log_num = log_comb(N - 1 - n_p, n - 1)
    log_den = log_comb(N - 1, n - 1)
    mu_p    = (n / N) * (1 - exp(log_num - log_den))
    
    logC_Nm2_ell = log_comb(N - 2, ell)
    logC_N_n     = log_comb(N,     n)
    log_common   = logC_Nm2_ell - logC_N_n
    common_fac   = exp(log_common)
    
    pi_mat = matrix(0, N, N)
    for (p in 1:(N - 1)) {
      Np_all = which(neigh[p, ])
      for (q in (p + 1):N) {
        if (neigh[p, q]) {
          pi = common_fac
        } else {
          Nq_all  = which(neigh[q, ])
          both    = intersect(Np_all, Nq_all)
          only_p  = setdiff(Np_all, c(q, both))
          only_q  = setdiff(Nq_all, c(p, both))
          
          a = length(only_p)
          b = length(only_q)
          w = length(both)
          c_cnt = N - 2 - a - b - w
          
          fail_p    = binom_prob(N - 2, ell, a + w)
          fail_q    = binom_prob(N - 2, ell, b + w)
          fail_both = binom_prob(N - 2, ell, a + b + w)
          
          prob_valid = 1 - fail_p - fail_q + fail_both
          
          pi = common_fac * prob_valid
        }
        pi_mat[p, q] = pi_mat[q, p] = pi
      }
    }
    var_G = (1 / n^2) * (
      sum(mu_p * (1 - mu_p)) +
        2 * sum(pi_mat[upper.tri(pi_mat)] -
                  outer(mu_p, mu_p, "*")[upper.tri(pi_mat)])
    )
    res[i] = pmax(var_G, 0)
  }
  res
}

exact_g = function(data, radii, mark_col, mark){
  out = exact_csr_rcpp2(data = data, radii = radii, mark_col = mark_col, mark = mark)
  out$`Empirical CSR Variance` = exact_G_variance(df = data, mark_col = mark_col, mark = mark, radii = radii)
  out = out %>%
    dplyr::mutate(z = (`Observed G` - `Empirical CSR`) / sqrt(`Empirical CSR Variance`),
                  z = ifelse(is.nan(z), 0, z),
                  p_value = 2 * (1 - pnorm(abs(z))))
  return(out)
}
