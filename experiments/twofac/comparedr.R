library(tidyverse)
library(brlavaan)
here::i_am("experiments/twofac/analysis.R")
load(here::here("experiments/simu_id.RData"))
load(here::here("experiments/simu_res_twofac.RData"))
simu_res <- simu_res_twofac

twofacpars <- c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")
mycols <- c(
  ML = "#E31A1C",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  Jackknife = "#B2DF8A",
  Bootstrap = "#33A02C",
  `Ozenne et al.` = "#FDBF6F",
  REML = "#FF7F00"
)

# ----- Prep results ----------------------------------------------------------

# Create (nested) results data frame, used for counting convergences etc.
res_nested <-
  simu_res |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows() |>
  drop_na() |>
  select(-starts_with("info"))

# This is the full (raw) results data frame
res <-
  res_nested |>
  mutate(param = lapply(truth, names), .before = est) |>
  unnest(param:truth) |>
  # # For the growth model, keep the first instance of param == "v"
  # filter(
  #   row_number() == 1,
  #   .by = c(simid, sim, dist, model, rel, n, method, param)
  # )
  # filter(param %in% twofacpars)
  mutate(type = case_when(
    grepl("=~", param) ~ "Lambda",
    grepl("fy~~fy|fx~~fx", param) ~ "Psi",
    grepl("~~", param) ~ "Theta",
    TRUE ~ "beta"
  )) |>
  mutate(
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

truth <-
  res |>
  summarise(
    across(truth, first),
    .by = c(model, param, rel)
  )

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
    seed,
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
  distinct(sim, seed, method, dist, n, rel, param, est, se, converged) |>
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
res_dr <- left_join(res_dr, select(simu_id, -seed))

# Rearrange columns
res_dr <-
  res_dr |>
  mutate(timing = NA, scaled_grad = NA, max_loglik = NA, optim_message = NA, Sigma_OK = NA, super_converged = NA) |>
  select(simid, seed, sim, dist, model, rel, n, method, param, est, se, truth,
         timing, converged, optim_message, super_converged) |>
  arrange(simid, sim, dist, model, rel, n, method, param) |>
  mutate(type = case_when(
    grepl("=~", param) ~ "Lambda",
    grepl("fy~~fy|fx~~fx", param) ~ "Psi",
    grepl("~~", param) ~ "Theta",
    TRUE ~ "beta"
  )) |>
  mutate(
    simid = as.integer(simid),
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

## ----- Tables ----------------------------------------------------------------
twofac_tab <-
  bind_rows(
    res,
    res_dr
  ) |>
  # mutate(converged = all(converged), .by = c(simid, seed)) |>
  # filter(converged) |>
  # filter(dist != "Kurtosis") |>
  # filter(param %in% twofacpars) |>
  summarise(
    mean_bias = mean(relbias, trim = 0.05),
    med_bias = median(relbias),
    rmse = sqrt(mean(relbias ^ 2, trim = 0.05)),
    coverage = mean(covered, na.rm = TRUE),
    .by = c(dist:method, type)
  ) |>
  select(-model)

twofac_tab |>
  arrange(dist, rel, n, type, method) |>
  print(n = 100)

## ----- Plots -----------------------------------------------------------------

p_biassamplesize_all <-
  res |>
  bind_rows(res_dr) |>
  # mutate(converged = all(converged), .by = c(simid, seed)) |>
  # filter(converged) |>
  filter(dist != "Kurtosis") |>
  # filter(param %in% twofacpars) |>
  summarise(
    bias = mean(relbias, trim = 0.05),
    .by = c(dist:method, type)
  ) |>
  mutate(
    n = as.numeric(factor(n)),
    rel = factor(rel, labels = c("Rel == 0.8", "Rel == 0.5"))
  ) |>
  ggplot(aes(n, bias, col = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggh4x::facet_nested(type ~ rel + dist, labeller = label_parsed) +
  scale_colour_manual(values = mycols) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  labs(x = "Sample size (n)", y = "Rel. mean bias", col = NULL) +
  guides(colour = guide_legend(nrow = 1, reverse = TRUE, position = "top")) +
  theme_bw(); p_biassamplesize_all

# save(
#   p_biassampsize_ours, p_conv_succ, p_dist_n50_all, p_perf_n50_all, p_perf_n50_type,
#   p1, p2, tab_conv, p_biassamplesize_all,
#   file = here::here("experiments/LATEST/simrestwofac.RData")
# )

bind_rows(
  res,
  res_dr
) |>
  filter(n == 50, rel == "Rel = 0.5") |>
  summarise(
    meanbias = mean(relbias, trim = 0.05),
    rmse = sqrt(mean(relbias ^ 2, trim = 0.05)),
    pu = mean(relbias < 0),
    coverage = mean(covered, na.rm = TRUE),
    .by = c(dist:method, type)
  ) |>
  mutate(var = rmse ^ 2 + meanbias ^ 2) |>
  pivot_longer(c(meanbias, rmse, pu, coverage), names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(
      metric,
      levels = c("meanbias", "var",  "rmse", "pu", "coverage"),
      labels = c("Mean", "Variance", "RMSE", "Prob. underest.", "Coverage")
    )
  ) |>
  ggplot(aes(value, type, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_vline(
    data = tibble(
      metric = factor(c("Mean", "RMSE", "Prob. underest.", "Coverage")),
      value = c(0, 0, 0.5, 0.95),
    ),
    aes(xintercept = value),
    linetype = "dashed"
  ) +
  scale_fill_manual(values = mycols) +
  ggh4x::facet_nested(rel + dist ~ metric, scales = "free", space = "free_y") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(labels = scales::percent),
      scale_x_continuous(labels = scales::percent),
      scale_x_continuous(limits = c(0.35, 0.65), labels = scales::percent),
      scale_x_continuous(limits = c(0.7, 1), labels = scales::percent)
    )
  ) +
  theme_bw() +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE, position = "bottom")) +
  labs(x = NULL, y = NULL, fill = NULL, title = "n = 50") +
  scale_y_discrete(labels = rev(c(
    expression(Theta),
    expression(Psi),
    expression(Lambda),
    expression(beta)
  )))
