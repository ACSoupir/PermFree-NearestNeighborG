# Regenerate the small synthetic example dataset shipped with the package.
#
# The proprietary mIF validation data used in the paper is deliberately NOT
# distributed with this package.  This synthetic pattern reproduces the small
# simulated example (30 points, 5 marked) used to demonstrate that the exact
# mean equals the average over all choose(30, 5) = 142506 subsets.
set.seed(20250611)
N <- 30L
n_marked <- 5L
sim_nnG <- data.frame(
  x = stats::runif(N),
  y = stats::runif(N),
  m = factor(c(rep("A", n_marked), rep("B", N - n_marked)))
)
stopifnot(nrow(sim_nnG) == 30L, sum(sim_nnG$m == "A") == 5L)
save(sim_nnG, file = "data/sim_nnG.rda", compress = "xz")
