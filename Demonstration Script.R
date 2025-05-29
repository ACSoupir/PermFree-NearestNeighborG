
# Study Script ------------------------------------------------------------


# Libraries ---------------------------------------------------------------

library(ggplot2)
library(magrittr)

source("functions.R")

# mIF Data ----------------------------------------------------------------
#10.1186/s13059-024-03435-z validation mIF
mif = readRDS("../../../../lab_fridley//Manley/Manley_SMI/data/UnivariateK_mif/manley_mif_validation_no-clinical.rds")
#get and clean example spatial image file
#ST RCC3_13 - clear cell - treatment naive - tumor
df = mif$spatial[[12]] %>%
  dplyr::mutate(x = (XMin + XMax)/2,
                y = (YMin + YMax)/2) %>%
  dplyr::select(x, y, 
                `FOXP3 (Opal 620) Positive Classification`,
                `CD8 (Opal 650) Positive Classification`, 
                `Classifier Label`, `Image Tag`) %>%
  dplyr::rename("FOXP3" = 3,
                "CD8" = 4) %>%
  dplyr::mutate(background = dplyr::case_when(rowSums(dplyr::select(., FOXP3, CD8)) > 1 |
                                                rowSums(dplyr::select(., FOXP3, CD8)) == 0 ~ 1,
                                              TRUE ~ 0),
                FOXP3 = ifelse(background == 1, 0, FOXP3),
                CD8 = ifelse(background == 1, 0, CD8)) %>%
  tidyr::pivot_longer(cols = c(FOXP3, CD8, background),
                      names_to = "label",
                      values_to = "Positive") %>%
  dplyr::filter(Positive == 1)

#plot to demonstrate 
mif_example = df %>%
  ggplot() + 
  geom_point(data = . %>%
               dplyr::filter(label == "background"),
             aes(x = x, y = y, shape = `Classifier Label`),
             size = 1) +
  geom_point(data = . %>%
               dplyr::filter(label != "background"),
             aes(x = x, y = y, fill = label),
             size = 1, color = "black", shape = 21) +
  scale_shape_manual(values = c(3, 16)) +
  scale_fill_manual(values= c("orange")) +
  guides(shape = guide_legend(title = "Compartment"),
         color = guide_legend(title = "")) +
  coord_equal() + 
  labs(title = unique(df[["Image Tag"]])) +
  facet_grid(~`Classifier Label`) +
  theme_classic() +
  theme(plot.title = element_text(hjust= 0.5)); mif_example
pdf("figures/example_tma.pdf", height = 5, width = 10)
mif_example
dev.off()

perm_example = permutations(data = df %>% dplyr::filter(`Classifier Label` == "Tumor"), 
             radii = 0:150, mark = "CD8", perms = 1000) %>%
  ggplot() +
  geom_line(aes(x = r, y = `Empirical CSR`, group = rep)) +
  geom_line(data = . %>%
              dplyr::group_by(r) %>%
              dplyr::summarise(t = mean(`Theoretical CSR`, na.rm = TRUE)),
            aes(x = r, y = t), color = "red") +
  geom_line(data = . %>%
              dplyr::group_by(r) %>%
              dplyr::summarise(t = mean(`Empirical CSR`, na.rm = TRUE)),
            aes(x = r, y = t), color = "green") +
  geom_line(data = . %>%
              dplyr::group_by(r) %>%
              dplyr::summarise(t = mean(`Observed G`, na.rm = TRUE)),
            aes(x = r, y = t), color = "orange") +
  labs(title = "Permutations vs Theoretical") +
  theme_classic() +
  theme(plot.title = element_text(hjust= 0.5)); perm_example
pdf("figures/permutations_example.pdf", height = 5, width = 10)
perm_example
dev.off()


# Compare methods ---------------------------------------------------------
#variables
r = 0:150
mark = "CD8"
dat = df %>%
  dplyr::filter(`Classifier Label` == "Tumor")
#functions
f_list = list(
  perm_res_10000 = function(){permutations(data = dat, radii = r, mark = mark, perms = 10000)},
  perm_res_1000 = function(){permutations(data = dat, radii = r, mark = mark, perms = 1000)},
  perm_res_100 = function(){permutations(data = dat, radii = r, mark = mark, perms = 100)},
  exact_res = function(){exact_csr(data = dat, radii = r, mark = mark)},
  exact_res_rcpp = function(){exact_csr_rcpp(data = dat, radii = r, mark = mark)},
  exact_res_rcpp2 = function(){exact_csr_rcpp2(data = dat, radii = r, mark = mark)}
)

#benchmarkign times
#100 samples per function
fun_names = rep(c("perm_res_10000", "perm_res_1000", "perm_res_100", "exact_res", "exact_res_rcpp", "exact_res_rcpp2"),
                100)
#randomize
set.seed(333)
fun_order = sample(seq(fun_names), size = length(fun_names), replace = FALSE)
fun_names = fun_names[fun_order]
#calculate all without preschedule to ensure no goofy things with order getting tied to cores
benchmark_dat = parallel::mclapply(fun_names, function(fn){
  microbenchmark::microbenchmark(f_list[[fn]](), times = 1) %>%
    data.frame()
}, mc.cores = 16, mc.preschedule = FALSE)
# saveRDS(benchmark_dat, "results_100iters.rds")
benchmark_dat = readRDS("results_100iters.rds")

#calculating actual values now
set.seed(333)
distr_dat = parallel::mclapply(f_list, function(f) f(), mc.cores = length(f_list))
# saveRDS(distr_dat, "results_distributions.rds")
distr_dat = readRDS("results_distributions.rds")

#exact benchmarks with memory
f_list_mem = list(
  perms10000 = function() { f_list$perm_res_10000(); return(NULL)},
  perms1000 = function() { f_list$perm_res_1000(); return(NULL)},
  perms100 = function() { f_list$perm_res_100(); return(NULL)},
  exact_csr = function() { f_list$exact_res(); return(NULL)},
  exact_csr_rcpp = function() { f_list$exact_res_rcpp(); return(NULL)},
  exact_csr_rcpp2 = function() { f_list$exact_res_rcpp2(); return(NULL)}
)

benchmark_mem = bench::mark(
  f_list_mem$perms10000(),
  f_list_mem$perms1000(),
  f_list_mem$perms100(),
  f_list_mem$exact_csr(),
  f_list_mem$exact_csr_rcpp(),
  f_list_mem$exact_csr_rcpp2(),
  iterations = 10
)
#saveRDS(benchmark_mem, "memory_benchmark.rds")
benchmark_mem = readRDS("memory_benchmark.rds")


# Plot Benchmarking times -------------------------------------------------

benchmark_dat2 = benchmark_dat %>%
  do.call(dplyr::bind_rows, .)
benchmark_dat2$fun = fun_names

computetimes = benchmark_dat2 %>%
  dplyr::mutate(time = time/1e9,
                Func = fun) %>%
  dplyr::mutate(Func = dplyr::case_when(Func == "perm_res_10000" ~ "Permuted CSR (10000)",
                                        Func == "perm_res_1000" ~ "Permuted CSR (1000)",
                                        Func == "perm_res_100" ~ "Permuted CSR (100)",
                                        Func == "exact_res" ~ "Exact CSR (R)",
                                        Func == "exact_res_rcpp" ~ "Exact CSR (R+Rcpp)",
                                        Func == "exact_res_rcpp2" ~ "Exact CSR (R+Rcpp Opt)"),
                Func = factor(Func, levels = c("Permuted CSR (100)", "Permuted CSR (1000)", "Permuted CSR (10000)",
                                               "Exact CSR (R)", "Exact CSR (R+Rcpp)", "Exact CSR (R+Rcpp Opt)"))) %>%
  ggplot() +
  geom_point(aes(x = Func, y = time), alpha = 0.2, size = 1, stroke=NA) +
  theme_classic() +
  labs(title = "Nearest Neighbor G Compute Time",
       x = "Functions",
       y = "Time (s)") +
  theme(plot.title = element_text(hjust= 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        text = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.text.x.bottom = element_text(color = "black")); computetimes
pdf("figures/nng_compute_time_comparison.pdf", height = 5, width = 7)
computetimes
dev.off()


# Plot curves -------------------------------------------------------------
X11()
distr_dat2 = distr_dat %>%
  dplyr::bind_rows(.id = "Func") %>%
  dplyr::mutate(Func = dplyr::case_when(Func == "perm_res_10000" ~ "Permuted CSR (10000)",
                                        Func == "perm_res_1000" ~ "Permuted CSR (1000)",
                                        Func == "perm_res_100" ~ "Permuted CSR (100)",
                                        Func == "exact_res" ~ "Exact CSR (R)",
                                        Func == "exact_res_rcpp" ~ "Exact CSR (R+Rcpp)",
                                        Func == "exact_res_rcpp2" ~ "Exact CSR (R+Rcpp Opt)"),
                Func = factor(Func, levels = c("Permuted CSR (100)", "Permuted CSR (1000)", "Permuted CSR (10000)",
                                               "Exact CSR (R)", "Exact CSR (R+Rcpp)", "Exact CSR (R+Rcpp Opt)")))
curves = distr_dat2 %>%
  dplyr::group_by(Func, r) %>%
  dplyr::select(-Group, -rep) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x))) %>%
  tidyr::pivot_longer(cols = c(`Theoretical CSR`, `Observed G`, `Empirical CSR`),
                      names_to = "Curve", values_to = "G(r)") %>%
  ggplot() + 
  geom_line(data = distr_dat$perm_res_1000,
            aes(x = r, y = `Empirical CSR`, group = rep), alpha = 0.05) +
  geom_line(aes(x = r, y = `G(r)`, color = Curve)) +
  scale_color_manual(values = c("Theoretical CSR" = "red", "Empirical CSR" = "green", "Observed G" = "orange")) +
  theme_classic() +
  facet_wrap(~Func) +
  guides(color = guide_legend(title = "")) +
  theme(legend.position = "bottom"); curves
pdf("figures/permuted_vs_exact.pdf", height = 6, width = 12)
curves
dev.off()


#plotting permutation 100 v others
tmp = distr_dat2 %>%
  dplyr::group_by(Func, r) %>%
  dplyr::select(-Group, -rep) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(.x))) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!grepl("Rcpp", Func))
distr_dat3 = dplyr::full_join(tmp %>%
                                dplyr::filter(Func == "Permuted CSR (100)") %>%
                                dplyr::select(r, `Empirical CSR`) %>%
                                dplyr::rename("100 Perms" = 2),
                              tmp %>%
                                dplyr::filter(Func != "Permuted CSR (100)") %>%
                                dplyr::select(Func, r, `Empirical CSR`)) %>%
  dplyr::mutate(`100 Perm Difference` = `Empirical CSR` - `100 Perms`)
curves2 = ggplot() + 
  geom_hline(yintercept = 0) +
  geom_line(data = distr_dat3 %>%
              dplyr::filter(grepl("Perm", Func)) %>% 
              dplyr::mutate(Func = as.character(Func)) %>%
              dplyr::rename("Permutations" = "Func"),
            aes(x = r, y = `100 Perm Difference`, color = Permutations)) +
  geom_line(data = distr_dat3 %>%
              dplyr::filter(!grepl("Perm", Func)),
            aes(x= r, y = `100 Perm Difference`, linetype = Func)) +
  facet_grid(Func ~. ) +
  theme_classic() +
  theme(legend.position = "bottom"); curves2
pdf("figures/permuted_vs_exact_sameplot.pdf", height = 5, width = 8)
curves2
dev.off()


# Make memory table -------------------------------------------------------
library(bench)
mem_pl = benchmark_mem %>%
  tidyr::unnest(c(time, gc)) %>% 
  dplyr::mutate(expression = gsub(".*\\$", "", as.character(expression))) %>%
  dplyr::mutate(expression = dplyr::case_when(expression == "perms10000()" ~ "Permuted CSR (10000)",
                                              expression == "perms1000()" ~ "Permuted CSR (1000)",
                                              expression == "perms100()" ~ "Permuted CSR (100)",
                                              expression == "exact_csr()" ~ "Exact CSR (R)",
                                              expression == "exact_csr_rcpp()" ~ "Exact CSR (R+Rcpp)",
                                              expression == "exact_csr_rcpp2()" ~ "Exact CSR (R+Rcpp Opt)"),
                expression = factor(expression, levels = c("Permuted CSR (100)", "Permuted CSR (1000)", "Permuted CSR (10000)",
                                               "Exact CSR (R)", "Exact CSR (R+Rcpp)", "Exact CSR (R+Rcpp Opt)"))) %>%
  ggplot() +
  geom_point(aes(x = mem_alloc, y = time, fill = expression), color = "black", shape = 21) +
  scale_color_bench_expr(scales::brewer_pal(type = "qual", palette = 3)) +
  theme_classic() +
  labs(title = "Memory and Time for Complete Spatial Randomness Calculations") +
  guides(fill = guide_legend(title = "Functions")) +
  theme(plot.title = element_text(hjust = 0.5)); mem_pl

pdf("figures/memory_usage.pdf", height = 6, width = 10)
mem_pl
dev.off()

#memory table
format_bytes <- function(bytes) {
  units <- c("B", "KB", "MB", "GB", "TB", "PB")
  idx <- pmin(floor(log(bytes, 1024)), length(units) - 1) 
  formatted_values <- sprintf("%.1f %s", bytes / (1024^idx), units[idx + 1])
  return(formatted_values)
}

mem_tab = benchmark_mem %>%
  tidyr::unnest(c(time, gc)) %>% 
  dplyr::mutate(expression = gsub(".*\\$", "", as.character(expression))) %>%
  dplyr::mutate(expression = dplyr::case_when(expression == "perms10000()" ~ "Permuted CSR (10000)",
                                              expression == "perms1000()" ~ "Permuted CSR (1000)",
                                              expression == "perms100()" ~ "Permuted CSR (100)",
                                              expression == "exact_csr()" ~ "Exact CSR (R)",
                                              expression == "exact_csr_rcpp()" ~ "Exact CSR (R+Rcpp)",
                                              expression == "exact_csr_rcpp2()" ~ "Exact CSR (R+Rcpp Opt)"),
                expression = factor(expression, levels = c("Permuted CSR (100)", "Permuted CSR (1000)", "Permuted CSR (10000)",
                                                           "Exact CSR (R)", "Exact CSR (R+Rcpp)", "Exact CSR (R+Rcpp Opt)"))) %>%
  dplyr::select(expression, time, mem_alloc) %>% 
  dplyr::group_by(expression) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(),
                                 list(mean = mean,
                                      sd = sd))) %>%
  dplyr::mutate(mem_alloc_mean = format_bytes(mem_alloc_mean),
                time_sd = paste0(round(time_sd * 100, 2), "ms")) %>%
  dplyr::rename("Function" = 1, "Time (mean)" = 2, "Time (sd)" = 3, "Memory (mean)" = 4, "Memory (sd)" = 5)
write.csv(mem_tab, "memory_table.csv")



# testing high permutation requirements -----------------------------------

permutation.mc = function(data, radii, mark, perms, cores){
  win = spatstat.geom::convexhull.xy(data[,1:2])
  pp = spatstat.geom::ppp(x = data$x,
                          y = data$y,
                          window = win,
                          marks = as.factor(data$label))
  G_obs = spatstat.explore::Gest(subset(pp, marks == mark), r = radii, rmax = max(radii), correction = "none") %>%
    as.data.frame() %>%
    dplyr::rename("Theoretical CSR" = 2, "Observed G" = 3)
  
  pcells = sum(data$label == mark)
  perm_list = lapply(split(1:perms, cut(seq_along(1:perms), cores, labels = FALSE)), length)
  perm = parallel::mclapply(seq_along(perm_list), function(perms2_n){
    lapply(1:perm_list[[perms2_n]], function(x){
      spatstat.explore::Gest(pp[sample(1:nrow(data), pcells, replace = F)], r =radii, rmax = max(radii), correction = "none") %>%
        as.data.frame() %>%
        dplyr::mutate(Group = "Permuted",
                      rep = c(0, unname(cumsum(do.call(c, perm_list))))[perms2_n] + x)
    }) %>%
      do.call(dplyr::bind_rows, .)%>%
      dplyr::rename("Theoretical CSR" = 2, "Empirical CSR" = 3)
  }, mc.cores = cores, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .)
  
  res = perm %>%
    dplyr::select(1,3,4,5) %>%
    # dplyr::group_by(r) %>%
    # dplyr::summarise(`Empirical CSR` = mean(raw, na.rm = TRUE)) %>%
    dplyr::full_join(G_obs, ., dplyr::join_by(r))
  return(res)
}

# set.seed(333)
# perms_100k = permutation.mc(data = dat, radii = r, mark = mark, perms = 1e5, cores = 16)
# perms_100k_summ = data.table::data.table(perms_100k)[,
#                                                      .(`Theoretical CSR` = mean(`Theoretical CSR`),
#                                                        `Observed G` = mean(`Observed G`),
#                                                        `Empirical CSR` = mean(`Empirical CSR`)), by = r]
# rm(perms_100k)
# saveRDS(perms_100k_summ, "perms_100k_summ.rds")
perms_100k_summ = readRDS("perms_100k_summ.rds")
# set.seed(333)
# perms_1m = permutation.mc(data = dat, radii = r, mark = mark, perms = 1e6, cores = 16)
# perms_1m_summ = data.table::data.table(perms_1m)[,
#                                                  .(`Theoretical CSR` = mean(`Theoretical CSR`),
#                                                    `Observed G` = mean(`Observed G`),
#                                                    `Empirical CSR` = mean(`Empirical CSR`)), by = r]
# rm(perms_1m)
# saveRDS(perms_1m_summ, "perms_1m_summ.rds")
perms_1m_summ = readRDS("perms_1m_summ.rds")

distr_dat4 = dplyr::full_join(tmp %>%
                                dplyr::filter(Func == "Permuted CSR (100)") %>%
                                dplyr::select(r, `Empirical CSR`) %>%
                                dplyr::rename("100 Perms" = 2),
                              tmp %>%
                                dplyr::filter(Func != "Permuted CSR (100)") %>%
                                dplyr::select(Func, r, `Empirical CSR`) %>%
                                dplyr::bind_rows(
                                  ., 
                                  dplyr::bind_rows(
                                    perms_100k_summ %>%
                                      dplyr::mutate(Func = "Permuted CSR (100k)") %>%
                                      dplyr::select(Func, r, `Empirical CSR`),
                                    perms_1m_summ %>%
                                      dplyr::mutate(Func = "Permuted CSR (1M)") %>%
                                      dplyr::select(Func, r, `Empirical CSR`)
                                  ))
                                ) %>%
  dplyr::mutate(`100 Perm Difference` = `Empirical CSR` - `100 Perms`)
curves3 = ggplot() + 
  geom_hline(yintercept = 0) +
  geom_line(data = distr_dat4 %>%
              dplyr::filter(grepl("Perm", Func)) %>% 
              dplyr::mutate(Func = as.character(Func)) %>%
              dplyr::rename("Permutations" = "Func"),
            aes(x = r, y = `100 Perm Difference`, color = Permutations)) +
  geom_line(data = distr_dat4 %>%
              dplyr::filter(!grepl("Perm", Func)),
            aes(x= r, y = `100 Perm Difference`, linetype = Func)) +
  facet_grid(Func ~. ) +
  theme_classic() +
  theme(legend.position = "bottom"); curves3
pdf("figures/permuted_vs_exact_sameplot_1Mperms.pdf", height = 5, width = 8)
curves3
dev.off()
