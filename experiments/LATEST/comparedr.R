# > colnames(res_filtered)
# [1] "simid"           "sim"             "dist"            "model"           "rel"
# [6] "n"               "method"          "param"           "est"             "se"
# [11] "truth"           "timing"          "converged"       "scaled_grad"     "max_loglik"
# [16] "Sigma_OK"        "optim_message"   "super_converged" "type"            "bias"
# [21] "relbias"         "covered"

dr_file1 <- here::here("experiments/GCM_est_combined_final.RData")
dr_file2 <- here::here("experiments/2FSEM_est_combined_final.RData")
if (!file.exists(dr_file1))
  download.file("https://osf.io/vjq5m/download", destfile = dr_file1)
if (!file.exists(dr_file2))
  download.file("https://osf.io/cw5b7/download", destfile = dr_file2)

# Growth model
# load("experiments/GCM_est_combined_final.RData")
# res_dr_growth <-
#   as_tibble(Results) |>
#   unite("simu", jobid, iteration, sep = ".") |>
#   mutate(
#     simu = as.integer(factor(simu)),
#     convergence = convergence == 1
#   ) |>
#   select(
#     sim = simu,
#     method,
#     dist,
#     n = nobs,
#     converged = convergence,
#     rel,
#     `est_Day0~~Day0`,
#     `est_i~~i`,
#     `est_i~1`,
#     `est_s~~s`,
#     `est_s~1`,
#     `est_i~~s`,
#     `se_Day0~~Day0`,
#     `se_i~~i`,
#     `se_i~1`,
#     `se_s~~s`,
#     `se_s~1`,
#     `se_i~~s`
#   ) |>
#   pivot_longer(
#     cols = c(starts_with("est_"), starts_with("se_")),
#     names_to = c(".value", "param"),
#     names_sep = "_"
#   ) |>
#   filter(method %in% c("JB", "BB", "Ozenne", "REML")) |>
#   distinct(sim, method, dist, n, rel, param, est, se, converged) |>
#   mutate(model = "growth")

# Two-factor SEM
load(here::here("experiments/2FSEM_est_combined_final.RData"))
res_dr_twofac <-
  as_tibble(Results) |>
  unite("simu", jobid, iteration, sep = ".") |>
  mutate(
    simu = as.integer(factor(simu)),
    convergence = convergence == 1
  ) |>
  select(
    sim = simu,
    method,
    dist,
    n = nobs,
    converged = convergence,
    rel,
    starts_with("est"),
    starts_with("se")
  ) |>
  select(-ends_with("~1")) |>
  pivot_longer(
    cols = c(starts_with("est_"), starts_with("se_")),
    names_to = c(".value", "param"),
    names_sep = "_"
  ) |>
  filter(method %in% c("JB", "BB", "Ozenne", "REML")) |>
  distinct(sim, method, dist, n, rel, param, est, se, converged) |>
  mutate(model = "twofac")

# Combine results
# res_dr <- bind_rows(res_dr_growth, res_dr_twofac)
res_dr <- res_dr_twofac

# Clean up
res_dr$dist <- factor(
  res_dr$dist,
  levels = c("Normal", "Kurtosis", "NonNormal"),
  labels = c("Normal", "Kurtosis", "Non-normal")
)
res_dr$dist <- as.character(res_dr$dist)

res_dr$rel <- factor(
  res_dr$rel,
  levels = c("REL80", "REL50"),
  labels = c(0.8, 0.5)
)
res_dr$rel <- as.numeric(as.character(res_dr$rel))

res_dr$n <- as.numeric(as.character(res_dr$n))

tmp_lev <- levels(res_dr$method)
tmp_lev <- gsub("JB", "Jackknife", tmp_lev)
tmp_lev <- gsub("BB", "Bootstrap", tmp_lev)
tmp_lev <- gsub("Ozenne", "Ozenne et al.", tmp_lev)
tmp_lev <- gsub("REML", "REML", tmp_lev)
levels(res_dr$method) <- tmp_lev

# Add truth and simid
levels(truth$rel) <- c(0.8, 0.5)
truth$rel <- as.numeric(as.character(truth$rel))
res_dr <- left_join(res_dr, truth)
res_dr <- left_join(res_dr, simu_id)

# Rearrange columns
res_dr <-
  res_dr |>
  mutate(timing = NA, scaled_grad = NA, max_loglik = NA, optim_message = NA, Sigma_OK = NA, super_converged = NA) |>
  select(simid, sim, dist, model, rel, n, method, param, est, se, truth, timing,
         converged, optim_message, super_converged) |>
  arrange(simid, sim, dist, model, rel, n, method, param) |>
  mutate(type = case_when(
    grepl("=~", param) ~ "Lambda",
    grepl("fy~~fy|fx~~fx", param) ~ "Psi",
    grepl("~~", param) ~ "Theta",
    TRUE ~ "beta"
  )) |>
  filter(sim <= 1000) |>
  mutate(
    simid = as.integer(simid),
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

## ----- Plots -----------------------------------------------------------------
p_biassamplesize_all <-
  bind_rows(
    res_filtered,
    res_dr
  ) |>
  filter(dist != "Kurtosis") |>
  summarise(
    bias = mean(relbias, trim = the_trim),
    .by = c(dist:method, type)
  ) |>
  mutate(
    n = as.numeric(factor(n)),
    rel = factor(rel, labels = c("Rel == 0.8", "Rel == 0.5"))
  ) |>
  ggplot(aes(n, bias, col = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggh4x::facet_nested(type ~ dist + rel, labeller = label_parsed) +
  scale_colour_manual(values = mycols) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample size (n)", y = "Rel. mean bias", col = NULL) +
  guides(colour = guide_legend(nrow = 1, reverse = TRUE, position = "top")) +
  theme_bw(); p_biassamplesize_all

save(
  p_biassampsize_ours, p_conv_succ, p_dist_n50_all, p_perf_n50_all, p_perf_n50_type,
  p1, p2, tab_conv, p_biassamplesize_all,
  file = here::here("experiments/LATEST/simrestwofac.RData")
)


# CHECK
DR_Res <-
  Results |>
  filter(dist == "NonNormal", rel == "REL50", method == "MLB") |>
  slice(1:5) |>
  select(rel, dist, method, nobs, seed, starts_with("est")) |>
  nest(est = starts_with("est")) |>
  mutate(est = map(est, \(x) {
    colnames(x) <- gsub("est_", "", colnames(x))
    x
  }))

i <- 1
# set.seed(DR_Res$seed[i])
dat <- gen_data_growth(n = 15, rel = 0.5, dist = "Non-normal", scale = 1, seed = DR_Res$seed[i])
mod <- txt_mod_growth(0.5)
fit <- growth(mod, dat, start = truth(dat), bounds = "standard", meanstructure = TRUE)
est <- coef(fit)
class(est) <- "numeric"
tinytest::expect_equal(est, unlist(DR_Res$est[i]))
