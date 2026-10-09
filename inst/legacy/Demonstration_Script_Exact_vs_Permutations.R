
# Simulation for exact permuted and exact CSR -----------------------------
library(magrittr)
library(ggplot2)
source("functions.R")
#simulating 30 points which is 17,100,720 permutations with 5 positive
set.seed(333)
sim_dat = data.frame(x = runif(30, min = 0, max = 10),
                     y = runif(30, min = 0, max = 10),
                     marks = sample(c(rep("pos", 5),
                                      rep("neg", 25)),
                                    30, replace = FALSE))
radii = seq(0, 15, 0.1)
mark = "pos"
cores = 16

#combination matrix
combo_matrix = combn(30, 5)
combo_df = as.data.frame(t(combo_matrix))

#running combos
sim_win = spatstat.geom::convexhull.xy(sim_dat[,c("x", "y")])
sim_pp = spatstat.geom::ppp(x = sim_dat$x,
                            y = sim_dat$y,
                            window = sim_win,
                            marks = sim_dat$marks)
sim_G_obs = spatstat.explore::Gest(subset(sim_pp, marks == mark), r = radii, rmax = max(radii), correction = "none") %>%
  as.data.frame() %>%
  dplyr::rename("Theoretical CSR" = 2, "Observed G" = 3)
sim_G_obs %>%
  ggplot() + geom_line(aes(x = r, y = `Observed G`))

pcells = sum(sim_dat$marks == mark)
perms = nrow(combo_df)
perm_list = split(1:perms, cut(seq_along(1:perms), cores, labels = FALSE))
sim_perm = parallel::mclapply(seq_along(perm_list), function(perms2_n){
  lapply(perm_list[[perms2_n]], function(x){
    spatstat.explore::Gest(sim_pp[as.numeric(combo_df[x,])], r =radii, rmax = max(radii), correction = "none") %>%
      as.data.frame() %>%
      dplyr::mutate(Group = "Permuted",
                    rep = x)
  }) %>%
    do.call(dplyr::bind_rows, .)%>%
    dplyr::rename("Theoretical CSR" = 2, "Empirical CSR" = 3)
}, mc.cores = cores, mc.preschedule = FALSE) %>%
  do.call(dplyr::bind_rows, .)

sim_res = sim_perm %>%
  dplyr::select(1,3,4,5) %>%
  dplyr::full_join(sim_G_obs, ., dplyr::join_by(r))



sim_summ_mean = data.table::data.table(sim_res)[,
                                           .(`Theoretical CSR` = mean(`Theoretical CSR`),
                                             `Observed G` = mean(`Observed G`),
                                             `Permuted CSR` = mean(`Empirical CSR`)), by = r]
exact_mean = exact_csr(data = sim_dat, mark_col = "marks",
                       radii = radii, mark = mark)


curves4 = sim_summ_mean %>%
  tidyr::pivot_longer(-r, names_to = "Metric", values_to = "G") %>%
  dplyr::mutate("Group" = "Permutation") %>%
  dplyr::bind_rows(., exact_mean %>%
                     dplyr::select(r, `Empirical CSR`, Group) %>%
                     dplyr::mutate(Metric = "Exact CSR") %>%
                     dplyr::rename("G" = "Empirical CSR"))  %>%
  dplyr::mutate(Metric = factor(Metric,
                                levels = c("Observed G", "Theoretical CSR",
                                           "Permuted CSR", "Exact CSR"))) %>%
  ggplot() + 
  geom_line(aes(x = r, y = G, color = Metric, linetype = Group), linewidth = 2) +
  theme_classic() +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_color_manual(values = c("lightblue", "lightgreen", "red", "purple")); curves4
# pdf("figures/all_permutations_vs_exact.pdf", height = 3, width = 8)
# curves4
# dev.off()

sim_summ_var = data.table::data.table(sim_res)[,
                                               .(`Theoretical CSR` = var(`Theoretical CSR`),
                                                 `Observed G` = var(`Observed G`),
                                                 `Permuted Variance` = var(`Empirical CSR`)), by = r]

exact_var = exact_G_variance(sim_dat, mark_col = "marks",
                             mark = "pos", 
                             r = radii)

dplyr::full_join(sim_summ_var,
                 data.frame(r = radii,
                            `Exact Variance` = exact_var, 
                            check.names = FALSE)) %>%
  dplyr::select(r, `Permuted Variance`, `Exact Variance`)  %>% View()


d = exact_g(data = sim_dat, radii = radii, mark_col = "marks", mark = mark)
d %>%
  ggplot() + 
  geom_line(aes(x = r, y = `Observed G`), color = 'black', linewidth = 2) +
  geom_line(aes(x = r, y = `Empirical CSR`), color = "blue") +
  geom_line(aes(x = r, y = `Empirical CSR` - sqrt(`Empirical CSR Variance`)), color = "red", alpha = 0.5) +
  geom_line(aes(x = r, y = `Empirical CSR` + sqrt(`Empirical CSR Variance`)), color = "red", alpha = 0.5) +
  geom_line(aes(x = r, y = `Observed G`, color = p_value), linewidth = 1)
